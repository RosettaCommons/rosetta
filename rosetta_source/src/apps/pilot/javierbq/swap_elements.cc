
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
public:

	NewMover():
		rebuild_max_iterations_(10),
		junction_start_(108),
		junction_end_(110),
		junction_ss_("LL"),
		junction_aa_("GG"),
		clash_dist_sq_( 3.5*3.5 )
	{	 }


	virtual void apply(Pose & pose) {
		Pose working_pose( pose );
		// Set up VarLenghtBuild to rebuild the junction
		SegmentRebuildOP seg_reb = new SegmentRebuild( 
										Interval( junction_start_, junction_end_ ),
										junction_ss_,
										junction_aa_,
										core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ),
										false //keep_known_bb_torsions_at_junctions
										);	
		BuildManager manager;
		manager.add( seg_reb );
		VarLengthBuild vlb( manager );

		// Rebuild the junction and let the transformation propagate through the structure
		Size iteration = 1;
		while( iteration <= rebuild_max_iterations_ ) {
			vlb.apply(working_pose);
			bool monomer_clash = check_monomer_bb_clash( working_pose );
			if( !monomer_clash ) {
				// If the backbone is not clashing align the new structure in the swap dimer configuration
				// check that remodeled region is not clashing with itself
					// make a new copy of the remodeled pose and align both copies to the input pose
					Pose new_pose( working_pose );	
					Pose ref_pose( pose );
					protocols::simple_moves::SuperimposeMover superimposer(working_pose);
					protocols::simple_moves::SuperimposeMover superimposer2(ref_pose);
					//superimposer.set_reference_pose( junction_end_, pose.total_residue() - 1);
					superimposer.set_target_range( junction_end_, pose.total_residue() - 1);
					superimposer.apply( ref_pose );
					ref_pose.dump_pdb("ref_pose.pdb");
					new_pose.dump_pdb("new_pose.pdb");
					working_pose.dump_pdb("working.pdb");
    			//superimposer2.set_reference_pose( 1, junction_start_ );
					superimposer2.set_target_range(1, junction_start_);
					superimposer2.apply( new_pose );
					ref_pose.dump_pdb("ref_pose_2.pdb");
					new_pose.dump_pdb("new_pose_2.pdb");
					working_pose.dump_pdb("working_2.pdb");
					// copy both copies into the input pose and return
					Pose dimer_pose(working_pose);
					// append the first residue with a jump
					dimer_pose.append_residue_by_jump( new_pose.residue(1), dimer_pose.total_residue() );
					// append the rest of the pose with bonds
					for(Size i = 2; i <= new_pose.total_residue(); i++) 
						dimer_pose.append_residue_by_bond( new_pose.residue(i) );
					bool dimer_clash = check_dimer_clash( dimer_pose );
					if( !dimer_clash ) {
						pose = dimer_pose;
						TR << "Swapped successfull after " << iteration << " iterations" << std::endl;
						set_last_move_status( protocols::moves::MS_SUCCESS );
						return;

					} else {
						TR << "Dimer is clashing after rebuild, trial " << iteration << " of " << rebuild_max_iterations_ << std::endl;
					}
			}	else {
				TR << "Monomer is clashing after rebuild, trial " << iteration << " of " << rebuild_max_iterations_ << std::endl;
			}
			iteration++;
		}
		TR << "Swap failed, maximum number of iterations tried" << std::endl;
		working_pose.dump_pdb("failed.pdb");
		set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
	}

	virtual std::string get_name() const {return "NewMover";}

private:
	bool check_monomer_bb_clash( Pose const & pose ) {
		for(Size i=1; i <= pose.total_residue(); ++i)
			for(Size j=i+1; j <= pose.total_residue(); ++j)
				if(pose.xyz(AtomID(pose.residue(i).atom_index("CA"),i)).distance_squared(pose.xyz(AtomID(pose.residue(j).atom_index("CA"),j))) < clash_dist_sq_) {
					TR << "monomer clash " << i << " " << j << std::endl;
					return true;
				}
		return false;
	}

	bool check_dimer_clash( Pose const & pose ) {
		const Size chain_length = pose.chain_sequence(1).size();
		for(Size i = 1; i <= chain_length ; ++i) 
			for(Size  j = 1 + chain_length; j <= pose.total_residue() ;  ++j) 
				if(pose.xyz(AtomID(pose.residue(i).atom_index("CA"),i)).distance_squared(pose.xyz(AtomID(pose.residue(j).atom_index("CA"),j))) < clash_dist_sq_)
					return true;
		return false;
	}

private:
	Size rebuild_max_iterations_;
	Size junction_start_;
	Size junction_end_;
	std::string junction_ss_;
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
