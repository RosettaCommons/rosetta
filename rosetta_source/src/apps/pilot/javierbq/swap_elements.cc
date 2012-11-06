
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   
/// @brief  
/// @author 


#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>

#include <numeric/xyz.functions.hh>
#include <basic/database/open.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/rms_util.hh>

#include <boost/lexical_cast.hpp>
#include <utility/exit.hh>
// Utility headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <protocols/simple_moves/symmetry/DetectSymmetryMover.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/SegmentRebuild.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/build/Interval.hh>
#include <core/chemical/ResidueType.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/picking_old/vall/util.hh>
#include <core/fragment/picking_old/FragmentLibraryManager.hh>
#include <core/fragment/Frame.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <numeric/xyzVector.hh>
#include <protocols/rigid/RB_geometry.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer TR("main");


class ThisApplication  {
public:
	ThisApplication();
	static void register_options();
private:
};

ThisApplication::ThisApplication()
{}

//OPT_1GRP_KEY( File, cluster, out )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
//  OPT(in::file::s);
}

class NewMover : public protocols::moves::Mover {
public:
	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pose::PoseOP PoseOP;
	typedef core::id::AtomID AtomID;
	typedef protocols::forge::build::BuildManager BuildManager;
	typedef protocols::forge::build::Interval Interval;
	typedef protocols::forge::build::SegmentRebuild SegmentRebuild;
	typedef protocols::forge::build::SegmentRebuildOP SegmentRebuildOP;
	typedef protocols::forge::components::VarLengthBuild VarLengthBuild;
	typedef core::kinematics::MoveMap MoveMap;
	typedef core::kinematics::MoveMapOP MoveMapOP;
	typedef core::fragment::FrameList FrameList;
	typedef core::fragment::FrameOP FrameOP;
	typedef core::fragment::FragData FragData;
	typedef numeric::xyzVector< Real > Vector;
public:

	NewMover():
		rebuild_max_iterations_(10),
		junction_start_(108),
		junction_end_(110),
		junction_ss_("LLLL"),
		junction_abego_("X"),
		junction_aa_("PPPP"),
		clash_dist_sq_( 3.5*3.5 ),
		nfrags_(100),
		frag_length_(5)
	{	 }

  /// @brief pick fragments of a given length, padding when necessary
  /// @param[in] complete_ss The complete secondary structure string, typically from a Pose.
  /// @param[in] complete_aa The complete amino acid string, typically from a Pose;
  ///            can be empty.  If empty, sequence bias is not used to pick fragments.
  /// @param[in] interval The interval [left, right] to pick fragments from; Pose
  ///  numbering (i.e. 1-based indexing).
  /// @param[in] frag_length The desired length of the fragments
  /// @param[in] n_frags The number of fragments to pick per position.
  FrameList pick_fragments(
  	std::string const & complete_ss,
  	std::string const & complete_aa,
  	utility::vector1< std::string > const & complete_abego,
  	Interval const & interval,
  	Size const frag_length,
  	Size const n_frags
  )
  {
  	using core::fragment::Frame;
  	using core::fragment::FrameOP;
  	using core::fragment::IndependentBBTorsionSRFD;
  
  	using core::fragment::picking_old::vall::pick_fragments;
  	using core::fragment::picking_old::vall::pick_fragments_by_ss;
  	using core::fragment::picking_old::vall::pick_fragments_by_ss_plus_aa;
  
  	FrameList frames;
  
  	for ( Size j = 0, je = interval.length(); j < je; ++j ) {
  		TR << "picking " << n_frags << " " << frag_length << "-mers for position " << ( interval.left + j ) << std::endl;
  
  		std::string ss_sub = complete_ss.substr( interval.left + j - 1, frag_length );
  		if ( ss_sub.length() < frag_length ) {
  			ss_sub.append( frag_length - ss_sub.length(), 'D' );
  		}
  
  		std::string aa_sub;
  		if ( !complete_aa.empty() ) {
  			aa_sub = complete_aa.substr( interval.left + j - 1, frag_length );
  			if ( aa_sub.length() < frag_length ) {
  				aa_sub.append( frag_length - aa_sub.length(), '.' );
  			}
  		} else {
  			aa_sub = "";
  		}
  		utility::vector1< std::string > abego_sub;
  		if ( complete_abego.size() > 0 ) {
  			runtime_assert( complete_ss.length() == complete_abego.size() );
  			Size pos( 1 );
  			abego_sub.resize( frag_length );
  			for( Size ii = interval.left + j; ii <= interval.left + j + frag_length - 1; ++ii, ++pos ) {
  				if ( ii > complete_abego.size() ) {
  					abego_sub[ pos ] = "X";
  				} else {
  					abego_sub[ pos ] = complete_abego[ ii ];
  				}
  			}
  		} else {
  			abego_sub.clear(); // make sure it is empty
  		}
  
  		FrameOP frame = new Frame( interval.left + j, frag_length );
  
  		frame->add_fragment( pick_fragments( ss_sub, aa_sub, abego_sub, n_frags, true, IndependentBBTorsionSRFD() ) );
  		frames.push_back( frame );
  	}
  
  	return frames;
  }

	bool clash_monomer(Pose const & pose) {
		for(Size i = 1; i <= pose.total_residue(); i++)
			for(Size j = i + 2; j <= pose.total_residue(); j++)
				if( pose.residue( i ).xyz("CA") - pose.residue( j ).xyz("CA") < 3.5)
					return true;
		return false;
	}

	virtual void apply(Pose & pose) {
		// pick up the fragments for the junction
		utility::vector1< std::string > abego;
		std::string ss = pose.secstruct();
		ss.replace(junction_start_ - 1, junction_end_ - junction_start_, junction_ss_);
		std::string sequence = pose.sequence();
		sequence.replace(junction_start_ - 1, junction_end_ - junction_start_, junction_aa_);
		Interval interval(junction_start_, junction_end_);
		FrameList frames = pick_fragments(ss, sequence, abego, interval, frag_length_, nfrags_);
		// make a pose with two copies of the input pose conected through a virutal atom
		Pose new_pose(pose);
		Vector com_A = protocols::geometry::center_of_mass(pose,1, junction_start_);
		Vector com_B = protocols::geometry::center_of_mass(pose,junction_end_ , pose.total_residue());

		// create the new residue
		core::chemical::ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );
		core::conformation::ResidueOP rsd( core::conformation::ResidueFactory::create_residue( rsd_set.name_map( "VRT" ) ) );
		rsd->set_xyz( "X", (com_A + com_B) / 2 );
		new_pose.append_residue_by_jump( *rsd, pose.total_residue() );

		new_pose.append_residue_by_jump( pose.residue(1), pose.total_residue(),"","",true );
		for(Size i=2; i <= pose.total_residue(); i++)
			new_pose.append_residue_by_bond( pose.residue(i) );
		// fix the backbone of half of one copy and half of the other.
		MoveMap movemap;
		movemap.set_jump(1,false);
		movemap.set_jump(2,false);
		for(Size i = junction_end_; i <= pose.total_residue() + junction_start_; i++) {
		//	movemap.set_bb(false);
			//movemap.set_chi(false);
		}
		
		TR << "inserting fragments" << std::endl;
		// connect a virtual atom from the center of mass of each moving region to the center of
		// mass of the other.
		// add a virtual atom in the center of mass of each rigid segment
		// make the fragment insertion in both junctions.

		bool clash = true;
		for(FrameList::iterator frame = frames.begin(); frame != frames.end(); frame++) {
			for(Size i = 1; i <= frames.size(); i++) {
				FrameOP frame = frames[i];
				for(Size j = 1; j <= frame->length(); j++) {
					FragData frag = frame->fragment(j);
					// insert fragment simultaneusly in the two junctions
					frag.apply(movemap, new_pose, junction_start_, junction_end_);
					frag.apply(movemap, new_pose, sequence.size() +  junction_start_ - 1, sequence.size() +  junction_end_ - 1);
					clash = false;// clash_monomer(new_pose);
					//Real rms = core::scoring::CA_rms(
					if( !clash ) 
						break;
					else
						TR << "monomer clash" << std::endl;
				}
				if( !clash ) break;
			}
			if( !clash ) break;
		}
		// check if the insertion move away the moving region and if it is close to the its partner
			//check for clashes
		pose = new_pose;
	}

	virtual std::string get_name() const { return "NewMover"; }

private:
	Size rebuild_max_iterations_;
	Size junction_start_;
	Size junction_end_;
	Size frag_length_;
	Size nfrags_;;
	std::string junction_ss_;
	std::string junction_abego_;
	std::string junction_aa_;
	Real clash_dist_sq_;
};

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

int main( int argc, char** argv ) {
	ThisApplication::register_options();
	devel::init( argc, argv );
	// mover
	protocols::moves::MoverOP protocol;
	protocol = new NewMover( );

	// run
	protocols::jd2::JobDistributor::get_instance()->go( protocol );
	
}
