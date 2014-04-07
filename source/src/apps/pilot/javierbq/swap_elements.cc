
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
#include <core/kinematics/FoldTree.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/picking_old/vall/util.hh>
#include <core/fragment/picking_old/FragmentLibraryManager.hh>
#include <core/fragment/Frame.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <fstream>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/util.hh>
#include <ObjexxFCL/format.hh>
#include <core/conformation/symmetry/SymmData.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

#include <utility/excn/Exceptions.hh>

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

// make 2 copies of the input pose, fix from one chain everything below the junction
// and in the ither chain everything above the junction and do fragment insertions in
// both sides. Then filter the structures where the moving elements are close to each other
class SwapElementsMover1 : public protocols::moves::Mover {
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

	SwapElementsMover1():
		rebuild_max_iterations_(50),
		junction_start_(108),
		junction_end_(110),
		junction_ss_("LLLL"),
		junction_abego_("X"),
		junction_aa_("PPPP"),
		clash_distance_( 3.5 ),
		distance_threshold_( 1.5 ),
		nfrags_(500),
		frag_length_(7)
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
			for(Size j =  1; j <= pose.total_residue(); j++) {
				if( i == j) continue;
				if( pose.residue( i ).xyz("CA") - pose.residue( j ).xyz("CA") < clash_distance_);
					return true;
			}
		return false;
	}

	virtual void apply(Pose & pose) {
		using protocols::moves::MS_SUCCESS;
		using protocols::moves::FAIL_DO_NOT_RETRY;
		// pick up the fragments for the junction
		utility::vector1< std::string > abego;
		std::string ss = pose.secstruct();
		ss.replace(junction_start_ - 1, junction_end_ - junction_start_, junction_ss_);
		std::string sequence = pose.sequence();
		sequence.replace(junction_start_ - 1, junction_end_ - junction_start_, junction_aa_);
		Interval interval(junction_start_, junction_end_);
		FrameList frames = pick_fragments(ss, sequence, abego, interval, frag_length_, nfrags_);
		// make a pose with two copies of the input pose conected through a virutal atom
		Pose working_pose(pose);
		Vector com_A = protocols::geometry::center_of_mass(pose,1, junction_start_);
		Size resi_nearest_to_com_A = protocols::geometry::return_nearest_residue( pose, 1, junction_start_, com_A);
		Vector com_B = protocols::geometry::center_of_mass(pose,junction_end_ , pose.total_residue());
		Size resi_nearest_to_com_B = protocols::geometry::return_nearest_residue( pose, junction_end_, pose.total_residue(), com_B);

		// create the new residue and locate it between the moving elements
		core::chemical::ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );
		core::chemical::ResidueType vrt( rsd_set.name_map( "VRT" ) ) ;
		core::conformation::ResidueOP anchor( core::conformation::ResidueFactory::create_residue( vrt ) ) ;
		anchor->set_xyz( "X", (com_A + com_B) / 2 );
		working_pose.append_residue_by_jump( *anchor, pose.total_residue() );

		working_pose.append_residue_by_jump( pose.residue(1), pose.total_residue(),"","",true );
		for(Size i=2; i <= pose.total_residue(); i++)
			working_pose.append_residue_by_bond( pose.residue(i) );
		// fix the backbone of half of one copy and half of the other.
		MoveMap movemap;
		movemap.set_jump(1,false);
		movemap.set_jump(2,false);
		for(Size i = 1; i <= working_pose.total_residue(); i++) {
			movemap.set_bb(true);
			movemap.set_chi(true);
		}

		TR << "inserting fragments" << std::endl;
		// connect a virtual atom from the center of mass of each moving region to the center of
		// mass of the other.
		core::conformation::ResidueOP vrt_a( core::conformation::ResidueFactory::create_residue( vrt ) );
		vrt_a->set_xyz( "X", com_B );
		working_pose.append_residue_by_jump( *vrt_a,  resi_nearest_to_com_A );
		Size idx_a = working_pose.fold_tree().downstream_jump_residue( 3 );
		core::conformation::ResidueOP vrt_b( core::conformation::ResidueFactory::create_residue( vrt ) );
		vrt_b->set_xyz( "X", com_A );
		working_pose.append_residue_by_jump( *vrt_b, pose.total_residue() +  resi_nearest_to_com_B );
		Size idx_b = working_pose.fold_tree().downstream_jump_residue( 4 );
		// add a virtual atom in the center of mass of each rigid segment
		// make the fragment insertion in both junctions.

		bool clash = true;
		Size trial = 1;
		for(FrameList::iterator frame = frames.begin(); frame != frames.end(); frame++) {
			for(Size i = 1; i <= frames.size(); i++) {
				FrameOP frame = frames[i];
				for(Size j = 1; j <= frame->length(); j++) {
					FragData frag = frame->fragment(j);
					// insert fragment simultaneusly in the two junctions
					frag.apply( movemap,  working_pose, junction_start_, junction_end_);
					frag.apply( movemap, working_pose, sequence.size() +  junction_start_ - 1, sequence.size() +  junction_end_ - 1);
					Real dist = working_pose.residue(idx_a).xyz("X").distance( working_pose.residue(idx_b).xyz("X") );
					TR << "distance " << dist << std::endl;
					clash =  clash_monomer(working_pose);
					if( !clash && dist < distance_threshold_ )
						break;
					else {
						TR << "monomer clash, trial "<< trial <<" of " << rebuild_max_iterations_ << std::endl;
						trial += 1;
						clash = ( trial <= rebuild_max_iterations_ ) ? true : false;
					}
				}
				if( !clash ) break;
			}
			if( !clash ) break;
		}
		// check if the insertion move away the moving region and if it is close to the its partner
			//check for clashes
		if( ! clash_monomer(working_pose) ) {
			set_last_move_status( MS_SUCCESS );
			pose = working_pose;
		} else {
			set_last_move_status( FAIL_DO_NOT_RETRY );
		}

	}

	virtual std::string get_name() const { return "SwapElementsMover1"; }

private:
	Size rebuild_max_iterations_;
	Size junction_start_;
	Size junction_end_;
	Size frag_length_;
	Size nfrags_;
	std::string junction_ss_;
	std::string junction_abego_;
	std::string junction_aa_;
	Real distance_threshold_;
	Real clash_distance_;
};

// insert fragments in one chain, perform alignment and rotation for
// getting the partner.
class SwapElementsMover2 : public protocols::moves::Mover {
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
	typedef numeric::xyzMatrix< Real > Matrix;
public:

	SwapElementsMover2():
		junction_start_(108),
		junction_end_(110),
		junction_ss_("LL"),
		junction_abego_("X"),
		junction_aa_("PP"),
		clash_distance_( 3.5 ),
		distance_threshold_( 3.0 ),
		nfrags_(5000),
		frag_length_(5),
		symm_min_( false )
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

	Size count_clashes(Pose const & pose) {
		Size clashes = 0;
		for(Size i = 1; i <= pose.total_residue(); i++) {
			for(Size j =  1; j <= pose.total_residue(); j++) {
				if( i == j) continue;
				if( pose.residue( i ).xyz("CA") - pose.residue( j ).xyz("CA") < clash_distance_)
					clashes++;
			}
		}
		return clashes;
	}

	virtual void apply(Pose & pose) {
		using protocols::moves::MS_SUCCESS;
		using protocols::moves::FAIL_DO_NOT_RETRY;
		using namespace ObjexxFCL::format;
		// pick up the fragments for the junction
		utility::vector1< std::string > abego;
		std::string ss = pose.secstruct();
		ss.replace(junction_start_ - 1, junction_end_ - junction_start_, junction_ss_);
		std::string sequence = pose.sequence();
		sequence.replace(junction_start_ - 1, junction_end_ - junction_start_, junction_aa_);
		Interval interval(junction_start_, junction_end_);
		FrameList frames = pick_fragments(ss, sequence, abego, interval, frag_length_, nfrags_);
		// make a copy of the pose
		Pose working_pose(pose);
		Vector com_A = protocols::geometry::center_of_mass(pose,1, junction_start_);
		Size resi_nearest_to_com_A = protocols::geometry::return_nearest_residue( pose, 1, junction_start_, com_A);
		Vector com_B = protocols::geometry::center_of_mass(pose,junction_end_ , pose.total_residue());
		Size resi_nearest_to_com_B = protocols::geometry::return_nearest_residue( pose, junction_end_, pose.total_residue(), com_B);

		// create the new residue and locate it between the moving elements
		core::chemical::ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );
		core::chemical::ResidueType vrt( rsd_set.name_map( "VRT" ) ) ;
		core::conformation::ResidueOP anchor( core::conformation::ResidueFactory::create_residue( vrt ) ) ;

		// connect a virtual atom from the center of mass of each segment and  to the center of
		// mass of each other.
		core::conformation::ResidueOP vrt_a( core::conformation::ResidueFactory::create_residue( vrt ) );
		vrt_a->set_xyz( "X", com_A );
		working_pose.append_residue_by_jump( *vrt_a,  resi_nearest_to_com_A );
		Size idx_a = working_pose.fold_tree().downstream_jump_residue( 1 );
		core::conformation::ResidueOP vrt_b( core::conformation::ResidueFactory::create_residue( vrt ) );
		vrt_b->set_xyz( "X", com_B );
		working_pose.append_residue_by_jump( *vrt_b, resi_nearest_to_com_B );
		Size idx_b = working_pose.fold_tree().downstream_jump_residue( 2 );

		core::conformation::ResidueOP center_to_vrt_a( core::conformation::ResidueFactory::create_residue( vrt ) );
		center_to_vrt_a->set_xyz( "X", (com_A + com_B)/2 );
		working_pose.append_residue_by_jump( *center_to_vrt_a, idx_a );
		Size idx_center_to_a = working_pose.fold_tree().downstream_jump_residue( 3 );

		core::conformation::ResidueOP center_to_vrt_b( core::conformation::ResidueFactory::create_residue( vrt ) );
		center_to_vrt_b->set_xyz( "X", (com_A + com_B)/2 );
		working_pose.append_residue_by_jump( *center_to_vrt_b, idx_b );
		Size idx_center_to_b = working_pose.fold_tree().downstream_jump_residue( 4 );


		// perform fragment insertion and align the vrt_a and vrt_b to the x axis and equidistant
		// to the y axis. After that make a copy of the pose, perform a 180 degrees rotation in the
		// z-axis and append the new pose
		Pose best_pose_after_fragment_insertion;
		Size best_pose_nclashes(999999);
		for(FrameList::iterator frame = frames.begin(); frame != frames.end(); frame++) {
			for(Size i = 1; i <= frames.size(); i++) {
				FrameOP frame = frames[i];
				for(Size j = 1; j <= frame->length(); j++) {
					Pose pose_for_fragment_insertion ( working_pose );
					FragData frag = frame->fragment(j);
					frag.apply( pose_for_fragment_insertion, junction_start_, junction_end_);
					// orient the modified pose
					const Vector vrt_a_pos = pose_for_fragment_insertion.residue( idx_a ).xyz("X");
					const Vector vrt_b_pos = pose_for_fragment_insertion.residue( idx_b ).xyz("X");
					const Vector center_vrt_a_vrt_b = ( vrt_a_pos + vrt_b_pos ) / 2;
					// new_center_to_a, new_center_to_b and center_vrt_a_vrt_b have to be aligned to the xy-plane,
					// in this way is is guarented that after a 180 degrees rotation around z-axis will generate the partner.

					// Align the center between the guides to the origin and align the guides to the xy-plane
					const Vector new_center_to_a_1 = pose_for_fragment_insertion.residue( idx_center_to_a ).xyz("X");
					const Vector new_center_to_b_1 = pose_for_fragment_insertion.residue( idx_center_to_b ).xyz("X");
					const Vector to_origin = ( new_center_to_a_1 + new_center_to_b_1 ) / 2;
  				Matrix id_rot_mat = numeric::xyzMatrix< core::Real >::identity();
					pose_for_fragment_insertion.apply_transform_Rx_plus_v( id_rot_mat, -1 *to_origin	);
					Real angle_vrt_a_orig_x_axis = numeric::angle_degrees(Vector(new_center_to_a_1[0],to_origin[1],new_center_to_a_1[2]), to_origin, Vector(to_origin[0] + 1, to_origin[1],to_origin[2]));
					Matrix y_rot = numeric::y_rotation_matrix_degrees(  angle_vrt_a_orig_x_axis );
					pose_for_fragment_insertion.apply_transform_Rx_plus_v( y_rot, Vector(0,0,0)	);
					// rotate around the z-axis and align the guides to the x-axis
					const Vector new_center_to_a_2 = pose_for_fragment_insertion.residue( idx_center_to_a ).xyz("X");
					const Vector new_center_to_b_2 = pose_for_fragment_insertion.residue( idx_center_to_b ).xyz("X");
					const Vector o2 = ( new_center_to_a_2 + new_center_to_b_2 ) / 2;
					Real angle_guide_a_origin_x_axis = numeric::angle_degrees(new_center_to_a_2, Vector(0,0,0), Vector(1,0,0));
					Matrix z_rot_x_align = numeric::z_rotation_matrix_degrees(  (new_center_to_a_2 < 0)? 1 : 1 * angle_guide_a_origin_x_axis );
					pose_for_fragment_insertion.apply_transform_Rx_plus_v( z_rot_x_align, Vector(0,0,0) );

					/*
					std::ofstream transformed( "transformed.pdb");
					core::io::pdb::dump_pdb(pose_for_rotation, transformed);
					const Vector p2a = pose_for_rotation.residue( idx_center_to_a ).xyz("X");
					const Vector p2b = pose_for_rotation.residue( idx_center_to_b ).xyz("X");
					const Vector c2a = pose_for_rotation.residue( idx_a ).xyz("X");
					const Vector c2b = pose_for_rotation.residue( idx_b ).xyz("X");
					const Vector center_2a_2b = (c2a + c2b) / 2;
					transformed << "HETATM  1    C   ACY     1    "<< F(8,3,p2a[0])  <<   F(8,3,p2a[1])<< F(8,3,p2a[2]) <<"  0.00  0.00" << std::endl;
					transformed << "HETATM  2    C   ACY     2    "<< F(8,3,p2b[0])  <<   F(8,3,p2b[1])<< F(8,3,p2b[2]) <<"  0.00  0.00" << std::endl;
					transformed << "HETATM  3    C   ACY     3    "<< F(8,3,center_2a_2b[0])  <<  F(8,3,center_2a_2b[1])<< F(8,3,center_2a_2b[2]) <<"  0.00  0.00" << std::endl;
					transformed << "HETATM  4    C   ACY     4    "<< F(8,3,c2a[0])  <<   F(8,3,c2a[1])<< F(8,3,c2a[2]) <<"  0.00  0.00" << std::endl;
					transformed << "HETATM  5    C   ACY     5    "<< F(8,3,c2b[0])  <<   F(8,3,c2b[1])<< F(8,3,c2b[2]) <<"  0.00  0.00" << std::endl;
					transformed.close();
					*/

					// rotate 180 degrees in the z axis to get the partner
					Pose pose_for_rotation = pose_for_fragment_insertion;
					Matrix z_rot = numeric::z_rotation_matrix_degrees( -180.0 );
					pose_for_rotation.apply_transform_Rx_plus_v( z_rot  , Vector(0,0,0)	);
					/*
					std::ofstream transformed_z( "transformed_z.pdb");
					core::io::pdb::dump_pdb(pose_for_rotation, transformed_z);
					const Vector p2a_z = pose_for_rotation.residue( idx_center_to_a ).xyz("X");
					const Vector p2b_z = pose_for_rotation.residue( idx_center_to_b ).xyz("X");
					const Vector c3a = pose_for_rotation.residue( idx_a ).xyz("X");
					const Vector c3b = pose_for_rotation.residue( idx_b ).xyz("X");
					const Vector center_3a_3b = (c3a + c3b) / 2;
					transformed_z << "HETATM  1    C   ACY     1   "<< F(8,3,p2a_z[0])  <<  F(8,3,p2a_z[1])<<  F(8,3,p2a_z[2]) <<"  0.00  0.00" << std::endl;
					transformed_z << "HETATM  2    C   ACY     2   "<< F(8,3,p2b_z[0])  <<  F(8,3,p2b_z[1])<<  F(8,3,p2b_z[2]) <<"  0.00  0.00" << std::endl;
					transformed_z << "HETATM  3    C   ACY     3   "<< F(8,3,center_3a_3b[0])  <<  F(8,3,center_3a_3b[1]) "<< F(8,3,center_3a_3b[3]) <<"  0.00  0.00" << std::endl;
					transformed_z << "HETATM  4    C   ACY     4   "<< F(8,3,c3a[0])  <<  F(8,3,c3a[1])<<  F(8,3,c3a[2]) <<"  0.00  0.00" << std::endl;
					transformed_z << "HETATM  5    C   ACY     5   "<< F(8,3,c3b[0])  <<  F(8,3,c3b[1])<<  F(8,3,c3b[2]) <<"  0.00  0.00" << std::endl;
					transformed_z.close();
					*/

    			// append pose_for_fragment_insertion and pose_for_rotation into the working pose
    //			core::pose::remove_virtual_residues( &pose_for_fragment_insertion );
    //			core::pose::remove_virtual_residues( &pose_for_rotation );
    			//core::pose::remove_upper_terminus_type_from_pose_residue( pose_for_rotation, 1);
    			//core::pose::remove_lower_terminus_type_from_pose_residue( two_chains_pose, two_chains_pose.total_residue() );
    			//two_chains_pose.append_residue_by_jump( pose_for_rotation.residue(1), two_chains_pose.total_residue(), "", "", true);

    			Pose two_chains_pose;
    			for(Size i = 1; i <= pose_for_fragment_insertion.total_residue(); i++) {
    				if( pose_for_fragment_insertion.residue( i ).name1() == 'X' ) continue;
    				two_chains_pose.append_residue_by_bond( pose_for_fragment_insertion.residue(i));
    			}
    			two_chains_pose.append_residue_by_jump( pose_for_rotation.residue(1), two_chains_pose.total_residue(), "", "", true);
    			for(Size i = 2; i <= pose_for_rotation.total_residue(); i++) {
    				if( pose_for_rotation.residue( i ).name1() == 'X' ) continue;
    				two_chains_pose.append_residue_by_bond( pose_for_rotation.residue(i));
    			}
					//eval the pose
					Size nclashes = count_clashes( two_chains_pose );
					if( nclashes < best_pose_nclashes ) {
						best_pose_after_fragment_insertion = two_chains_pose;
						best_pose_nclashes = nclashes;
					}

				}
			}
		}

		TR << "best dimer number of clashes " << best_pose_nclashes << std::endl;
		// Set the movemap, add constraints and perform symmetric minimization to the pose
		if( symm_min_ ) {
			TR << "Performing a symmetric minimization of the swapped dimer" << std::endl;
			// make a symmetric version of the pose
			Pose symm_pose = best_pose_after_fragment_insertion.split_by_chain( 1 );
			const Size mono_size = symm_pose.total_residue();
			// set up the movemap
			core::scoring::ScoreFunctionCOP scorefxn = new core::scoring::symmetry::SymmetricScoreFunction( *core::scoring::getScoreFunction() );
			const Size junction_new_end = junction_end_ - (junction_aa_.size() - (junction_end_ - junction_start_ + 1) );
			// set the junction as movable
			MoveMapOP movemap = new MoveMap;
			for(Size i = junction_start_; i <= junction_new_end; i++) {
				movemap->set_bb(i, true);
				movemap->set_chi(i,true);
			}
			/*
			Vector com_A_chA = protocols::geometry::center_of_mass(symm_pose,1, junction_start_);
			Size resi_nearest_to_com_A_chA = protocols::geometry::return_nearest_residue( symm_pose, 1, junction_start_, com_A_chA);
			Vector com_B_chA = protocols::geometry::center_of_mass(symm_pose,junction_new_end, mono_size);
			Size resi_nearest_to_com_B_chA = protocols::geometry::return_nearest_residue( symm_pose, 1, junction_start_, com_B_chA);

			Vector com_A_chB = protocols::geometry::center_of_mass(symm_pose, mono_size + 1, mono_size + junction_start_);
			Size resi_nearest_to_com_A_chB = protocols::geometry::return_nearest_residue( symm_pose, mono_size + 1, mono_size + junction_start_, com_A_chB);
			Vector com_B_chB = protocols::geometry::center_of_mass(symm_pose,mono_size + junction_new_end, 2*mono_size );
			Size resi_nearest_to_com_B_chB = protocols::geometry::return_nearest_residue( symm_pose, mono_size + junction_new_end, 2*mono_size, com_B_chA);
			 //add guides and constraints between them
			core::conformation::ResidueOP vrt_A_chA( core::conformation::ResidueFactory::create_residue( vrt ) );
			core::conformation::ResidueOP vrt_A_chB( core::conformation::ResidueFactory::create_residue( vrt ) );
			core::conformation::ResidueOP vrt_B_chA( core::conformation::ResidueFactory::create_residue( vrt ) );
			core::conformation::ResidueOP vrt_B_chB( core::conformation::ResidueFactory::create_residue( vrt ) );
			vrt_A_chA->set_xyz( "X", com_A_ch_A );
			vrt_A_chB->set_xyz( "X", com_A_ch_B );
			vrt_B_chA->set_xyz( "X", com_B_ch_A );
			vrt_B_chB->set_xyz( "X", com_B_ch_B );
			working_pose.append_residue_by_jump( *vrt_A_chA,  resi_nearest_to_com_A_chA ); //jump 2?
			working_pose.append_residue_by_jump( *vrt_A_chB,  resi_nearest_to_com_A_chB );
			working_pose.append_residue_by_jump( *vrt_B_chA,  resi_nearest_to_com_B_chA );
			working_pose.append_residue_by_jump( *vrt_B_chB,  resi_nearest_to_com_B_chB );
			// check that the jumps are write create the constriants and add them to symm_pose
			//Size idx_a = working_pose.fold_tree().downstream_jump_residue( 1 );


			std::string db_file = "symmetry/cyclic/C2_Z.sym";
			std::string path_to_symdef = basic::database::full_name(db_file);
			core::conformation::symmetry::SymmData symmdef;
			symmdef.read_symmetry_data_from_file( path_to_symdef );
			core::pose::symmetry::make_symmetric_pose( symm_pose, symmdef );



			//minimize
			option[OptionKeys::symmetry::symmetry_definition].value( path_to_symdef );
			core::pose::symmetry::make_symmetric_movemap( symm_pose, *movemap );
			protocols::simple_moves::symmetry::SymMinMover symm_min( movemap, scorefxn, "dfpmin_armijo_nonmonotone", 1e-5, true, false, false );
			symm_min.apply(symm_pose);
			option[OptionKeys::symmetry::symmetry_definition].value( "" );

			pose = core::pose::symmetry::get_asymmetric_pose_copy_from_symmetric_pose( symm_pose );
*/

		} else {
			pose = best_pose_after_fragment_insertion;
		}

	//	if( ! clash_monomer(pose_for_fragment_insertion) ) {
	//		set_last_move_status( MS_SUCCESS );
//			pose = working_pose;
	//	} else {
	//		set_last_move_status( FAIL_DO_NOT_RETRY );
	//	}


	}

	virtual std::string get_name() const { return "SwapElementsMover2"; }

private:
	Size junction_start_;
	Size junction_end_;
	Size frag_length_;
	Size nfrags_;
	std::string junction_ss_;
	std::string junction_abego_;
	std::string junction_aa_;
	Real distance_threshold_;
	Real clash_distance_;
	bool symm_min_;
};

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

int main( int argc, char** argv ) {
	try {

	ThisApplication::register_options();
	devel::init( argc, argv );
	// mover
	protocols::moves::MoverOP protocol;
	protocol = new SwapElementsMover2( );

	// run
	protocols::jd2::JobDistributor::get_instance()->go( protocol );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
