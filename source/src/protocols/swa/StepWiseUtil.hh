// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Util.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_StepWiseUtil_HH
#define INCLUDED_protocols_swa_StepWiseUtil_HH


#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <string>
#include <map>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>

//Auto Headers
#include <core/id/AtomID.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>



using namespace core;
typedef  numeric::xyzMatrix< Real > Matrix;

namespace protocols {
namespace swa {

	typedef std::map< std::string, pose::PoseOP > PoseList;

	Size
	make_cut_at_moving_suite( pose::Pose & pose, Size const & moving_suite );

	bool
	is_cutpoint_closed( pose::Pose const & pose, Size const seq_num );

	utility::vector1< Size >
	get_cutpoint_closed( pose::Pose const & pose );

	bool
	Is_close_chain_break(pose::Pose const & pose);

	// Undefined, commenting out to fix PyRosetta build  bool Contain_seq_num(Size const & seq_num, utility::vector1< Size > const & residue_list);

	//following copies code that is in rna/StepWiseRNA_Util? Remove the latter?
	void
	Output_boolean(std::string const & tag, bool boolean, std::ostream & TR );

	void
	Output_boolean(bool boolean, std::ostream & TR );

	void
	Output_movemap(kinematics::MoveMap const & mm, Size const total_residue, std::ostream & TR);

	void
	Figure_out_moving_residues( kinematics::MoveMap & mm, pose::Pose const & pose,
															utility::vector1< Size > const & fixed_res,
															bool const move_takeoff_torsions = true,
															bool const move_jumps_between_chains = false
															);

	void
	pdbslice( core::pose::Pose & new_pose,
						core::pose::Pose const & pose,
						utility::vector1< core::Size > const & slice_res );

	void
	pdbslice( pose::Pose & pose,
						utility::vector1< Size > const & slice_res );


	/// @brief  Superimpose mod_pose onto ref_pose using the mapping of residues from
	/// mod_pose to ref_pose given by res_map
	Real
	superimpose_pose(
									 pose::Pose & mod_pose,
									 pose::Pose const & ref_pose,
									 std::map< Size, Size > const & res_map
									 );

	id::AtomID_Map< id::AtomID >
 	create_alignment_id_map(	pose::Pose const & mod_pose, pose::Pose const & ref_pose,
													utility::vector1< Size > const & superimpose_res );

 	id::AtomID_Map< id::AtomID >
 	create_alignment_id_map(	pose::Pose const & mod_pose,
													pose::Pose const & ref_pose,
													std::map< core::Size, core::Size > res_map );

	utility::vector1< Size > const
	convert_to_working_res( utility::vector1< Size > const & res_vector,
													utility::vector1< Size > const & working_res ) ;


	scoring::constraints::ConstraintSetOP
	constraint_set_slice( scoring::constraints::ConstraintSetOP & cst_set,
												utility::vector1< Size > const & slice_res,
												pose::Pose const & pose,
												pose::Pose const & full_pose );

	utility::vector1< std::string > load_s_and_l();


	std::string
	get_file_name( std::string const & silent_file, std::string const & tag );

	void
	check_scorefxn_has_constraint_terms_if_pose_has_constraints( pose::Pose const & pose,
																															 scoring::ScoreFunctionOP & scorefxn );

	utility::vector1< Size >
	merge_vectors( utility::vector1< Size > const & vec1,
								 utility::vector1< Size > const & vec2 );

	void
	get_euler_angles( Real & alpha, Real & beta, Real & gamma, Matrix M1, Matrix M2, bool const verbose = true );

	void
	create_euler_rotation(
												Matrix & M,
												Real const & alpha,
												Real const & beta,
												Real const & gamma,
												Vector const & /* axis1 not actually used*/,
												Vector const & axis2,
												Vector const & axis3
												);

	void
	create_euler_rotation(
												Matrix & M,
												Real const & alpha,
												Real const & beta,
												Real const & gamma );

	void
	translate( pose::Pose & pose, Vector const shift,
						 pose::Pose const & ref_pose,
						 utility::vector1< Size > const & moving_res );

	void
	rotate( pose::Pose & pose, Matrix const M,
					pose::Pose const & ref_pose,
					utility::vector1< Size > const & moving_res,
					Vector const & centroid );

	void
	rotate( pose::Pose & pose, Matrix const M,
					pose::Pose const & ref_pose,
					utility::vector1< Size > const & moving_res );

	void
	get_base_centroid_and_rotation_matrix( pose::Pose const & pose, Size const i, Vector & centroid, Matrix & M );

	void
	translate_and_rotate_residue_to_origin( pose::Pose & pose, Size const i, utility::vector1< Size > const moving_res, bool const do_not_rotate = false );

	void
	translate_and_rotate_residue_to_origin( pose::Pose & pose, Size const i );

	void
	add_to_atom_id_map_after_checks( std::map< id::AtomID, id::AtomID> & atom_id_map,
																	 std::string const & atom_name,
																	 Size const & n1, Size const & n2,
																	 pose::Pose const & pose1, pose::Pose const & pose2 );

	Real
	superimpose_at_fixed_res_and_get_all_atom_rmsd( pose::Pose & pose, pose::Pose const & native_pose );

	Size
	get_number_missing_residue_connections( pose::Pose & pose );

	bool
	is_at_terminus( core::pose::Pose const & pose, Size const i );

	void
	merge_in_other_pose( pose::Pose & pose, pose::Pose const & pose2, Size const & merge_res );

	utility::vector1< Size >
	merge_two_poses_using_full_model_info( pose::Pose & pose,
																				 pose::Pose const & pose1,
																				 pose::Pose const & pose2,
																				 Size const & merge_res );

	utility::vector1< Size >
	merge_two_poses( pose::Pose & pose,
									 pose::Pose const & pose1,
									 pose::Pose const & pose2,
									 utility::vector1< Size > const & working_res1,
									 utility::vector1< Size > const & working_res2,
									 Size const & merge_res,
									 bool const fix_first_pose = true );

	void
	slice( pose::Pose & sliced_out_pose,
				 pose::Pose const & pose,
				 utility::vector1< Size > const & slice_res );

	void
	slice_out_pose( pose::Pose & pose,
									pose::Pose & sliced_out_pose,
									utility::vector1< Size > const & residues_to_delete );

	void
	correctly_add_cutpoint_variants( core::pose::Pose & pose, Size const res_to_add, 	bool const check_fold_tree = true );


	void
	fix_up_residue_type_variants_at_strand_beginning( pose::Pose & pose, Size const res );

	void
	fix_up_residue_type_variants_at_strand_end( pose::Pose & pose, Size const res );

	void
	fix_up_residue_type_variants( pose::Pose & pose );

	void
	switch_focus_to_other_pose( pose::Pose & pose, Size const & focus_pose_idx );

	core::pose::PoseOP
	get_pdb_and_cleanup( std::string const input_file,
											 core::chemical::ResidueTypeSetCAP rsd_set );

	void
	cleanup( pose::Pose & pose );

	void
	get_other_poses( 	utility::vector1< pose::PoseOP > & other_poses,
										utility::vector1< std::string > const & other_files,
										core::chemical::ResidueTypeSetCAP rsd_set );

	utility::vector1< Size >
	figure_out_moving_chain_break_res( pose::Pose const & pose, kinematics::MoveMap const & mm );


}
}

#endif
