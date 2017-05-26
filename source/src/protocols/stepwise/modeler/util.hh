// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/stepwise/modeler/util.hh
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_StepWiseUtil_HH
#define INCLUDED_protocols_stepwise_StepWiseUtil_HH

#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/full_model_info/FullModelInfo.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.fwd.hh>

//Auto Headers
#include <core/id/AtomID.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

/////////////////////////////////////////////////////////////////////////////////
// All functions for cutting/pasting/merging poses that happens
//  in stepwise monte carlo or assembly.
// Try to put other functions in output_util.hh, scoring_util.hh, etc.
/////////////////////////////////////////////////////////////////////////////////


namespace protocols {
namespace stepwise {
namespace modeler {

typedef std::map< std::string, core::pose::PoseOP > PoseList;

core::Size
make_cut_at_moving_suite( core::pose::Pose & pose, core::Size const & moving_suite );

core::Size
make_cut_at_moving_suite( core::kinematics::FoldTree & fold_tree, core::Size const & moving_suite );

core::Size
find_jump_number_at_suite( core::kinematics::FoldTree const & fold_tree, core::Size const & moving_suite );

core::Size
look_for_unique_jump_to_moving_res( core::kinematics::FoldTree const & fold_tree, core::Size const & i );

bool
is_cutpoint_closed( core::pose::Pose const & pose, core::Size const seq_num );

utility::vector1< core::Size >
get_cutpoint_closed( core::pose::Pose const & pose );

utility::vector1< core::Size >
merge_vectors( utility::vector1< core::Size > const & vec1,
	utility::vector1< core::Size > const & vec2 );

void
merge_in_other_pose_by_bond( core::pose::Pose & pose, core::pose::Pose const & pose2, core::Size const merge_res );

void
merge_in_other_pose_by_jump( core::pose::Pose & pose, core::pose::Pose const & pose2,
	core::Size const lower_merge_res, core::Size const upper_merge_res );

void
merge_in_other_pose( core::pose::Pose & pose, core::pose::Pose const & pose2,
	core::Size const lower_merge_res, core::Size const upper_merge_res,
	bool const connect_residues_by_bond );

utility::vector1< core::Size >
merge_two_poses_using_full_model_info( core::pose::Pose & pose,
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	core::Size const lower_merge_res,
	core::Size const upper_merge_res,
	bool const connect_residues_by_bond );

utility::vector1< core::Size >
merge_two_poses( core::pose::Pose & pose,
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	utility::vector1< core::Size > const & working_res1,
	utility::vector1< core::Size > const & working_res2,
	core::Size const lower_merge_res,
	core::Size const upper_merge_res,
	bool const connect_residues_by_bond,
	bool const fix_first_pose = true );

void
declare_chemical_bonds_at_cutpoints( core::pose::Pose & pose,
	core::pose::Pose const & source_pose,
	utility::vector1< core::Size > const & working_res,
	utility::vector1< core::Size > const & source_working_res );

void
declare_chemical_bonds_at_cutpoints( core::pose::Pose & pose,
	core::pose::Pose const & source_pose,
	utility::vector1< core::Size > const & working_res );

void
slice( core::pose::Pose & sliced_out_pose,
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & slice_res );

void
slice_out_pose( core::pose::Pose & pose,
	core::pose::Pose & sliced_out_pose,
	utility::vector1< core::Size > const & residues_to_delete );

core::Size
check_jump_to_previous_residue_in_chain( core::pose::Pose const & pose, core::Size const i,
	utility::vector1< core::Size > const & current_element );

core::Size
check_jump_to_previous_residue_in_chain( core::pose::Pose const & pose, core::Size const i,
	utility::vector1< core::Size > const & current_element,
	core::pose::full_model_info::FullModelInfo const & full_model_info );

core::Size
check_jump_to_previous_residue_in_chain( core::pose::Pose const & pose, core::Size const i,
	utility::vector1< core::Size > const & current_element,
	utility::vector1< core::Size > const & res_list,
	utility::vector1< core::Size > const & chains_in_full_model );

core::Size
check_jump_to_next_residue_in_chain( core::pose::Pose const & pose, core::Size const i,
	utility::vector1< core::Size > const & current_element );

core::Size
check_jump_to_next_residue_in_chain( core::pose::Pose const & pose, core::Size const i,
	utility::vector1< core::Size > const & current_element,
	core::pose::full_model_info::FullModelInfo const & full_model_info );

core::Size
check_jump_to_next_residue_in_chain( core::pose::Pose const & pose, core::Size const i,
	utility::vector1< core::Size > const & current_element,
	utility::vector1< core::Size > const & res_list,
	utility::vector1< core::Size > const & chains_in_full_model );

void
switch_focus_to_other_pose( core::pose::Pose & pose,
	core::Size const & focus_pose_idx,
	core::scoring::ScoreFunctionCOP scorefxn = 0 );

bool
switch_focus_among_poses_randomly( core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP scorefxn = 0,
	bool force_switch = false );

utility::vector1< core::Size >
figure_out_moving_chain_break_res( core::pose::Pose const & pose, core::kinematics::MoveMap const & mm );

bool
check_for_input_domain( core::pose::Pose const & pose,
	utility::vector1< core::Size > const & partition_res );

bool
check_for_input_domain( core::pose::Pose const & pose );

void make_variants_match(
	core::pose::Pose & pose,
	core::pose::Pose const & reference_pose,
	core::Size const n,
	core::chemical::VariantType const variant_type );

utility::vector1< core::Size >
figure_out_moving_cutpoints_closed( core::pose::Pose const & pose,
	utility::vector1< core::Size > const & moving_partition_res );

utility::vector1< core::Size >
figure_out_moving_cutpoints_closed_from_moving_res( core::pose::Pose const & pose, core::Size const moving_res );

utility::vector1< core::Size >
figure_out_moving_cutpoints_closed_from_moving_res( core::pose::Pose const & pose, utility::vector1< core::Size > const & moving_res_list );

void
figure_out_moving_chain_breaks( core::pose::Pose const & pose,
	utility::vector1< core::Size > const & moving_partition_res,
	utility::vector1< core::Size > & cutpoints_closed,
	utility::vector1< core::Size > & five_prime_chain_breaks,
	utility::vector1< core::Size > & three_prime_chain_breaks,
	utility::vector1< core::Size > & chain_break_gap_sizes );

core::Size
figure_out_reference_res_for_suite( core::pose::Pose const & pose, core::Size const moving_res );

utility::vector1< bool >
get_partition_definition( core::pose::Pose const & pose, core::Size const & moving_suite );

utility::vector1< bool >
get_partition_definition_by_jump( core::pose::Pose const & pose, core::Size const & jump_nr /*jump_number*/ );

void
reroot_based_on_full_model_info( core::pose::Pose & pose );

void
reroot_based_on_full_model_info( core::pose::Pose & pose,
	utility::vector1< core::Size > const & root_partition_res );

utility::vector1< core::Size >
figure_out_moving_partition_res_for_suite( core::pose::Pose const & pose,
	core::Size const moving_res,
	core::Size const reference_res );

utility::vector1< core::Size >
figure_out_moving_partition_res_for_jump( core::pose::Pose const & pose,
	core::Size const jump_nr );

void
figure_out_root_and_moving_partition_res( core::pose::Pose const & pose, core::Size const moving_res,
	utility::vector1< core::Size > & root_partition_res,
	utility::vector1< core::Size > & moving_partition_res );

utility::vector1< core::Size >
figure_out_moving_partition_res( core::pose::Pose const & pose,
	core::Size const moving_res );

utility::vector1< core::Size >
figure_out_moving_partition_res( core::pose::Pose const & pose,
	utility::vector1< core::Size > const & moving_res_list );

utility::vector1< core::Size >
figure_out_root_partition_res( core::pose::Pose const & pose,
	utility::vector1< core::Size > const & moving_res_list );

bool
revise_root_and_moving_res( core::pose::Pose & pose, core::Size & moving_res /* note that this can change too*/ );

bool
revise_root_and_moving_res_list( core::pose::Pose & pose,
	utility::vector1< core::Size > & moving_res_list /* note that this can change too*/ );

core::Size
find_downstream_connection_res( core::pose::Pose const & pose,
	utility::vector1< core::Size > const & moving_partition_res );

core::Size
split_pose( core::pose::Pose & pose, core::Size const moving_res, core::Size const reference_res );

void
split_pose( core::pose::Pose & pose, utility::vector1< core::Size > const & moving_res_list );

void
fix_up_jump_atoms( core::pose::Pose & pose );

void
fix_up_jump_atoms_and_residue_type_variants( core::pose::Pose & pose_to_fix );

void
fix_protein_jump_atom( core::pose::Pose & pose, core::Size const res, std::string const & atom_name );

void
add_to_pose_list( utility::vector1< core::pose::PoseOP > & pose_list, core::pose::Pose const & pose, std::string const & pose_tag );

// deprecate soon -- after protein/RNA unification
bool
is_protein( core::pose::Pose const & pose, utility::vector1< core::Size > const & moving_res_list );

utility::vector1< core::Size >
get_domain_boundary_suites( core::pose::Pose const & pose );

utility::vector1< core::Size >
get_domain_boundary_res( core::pose::Pose const & pose );

utility::vector1< core::Size >
get_moving_res_including_domain_boundaries( core::pose::Pose const & pose, utility::vector1< core::Size > const & moving_res_list );

utility::vector1< core::Size >
get_all_working_moving_res( working_parameters::StepWiseWorkingParametersCOP working_parameters);

void
virtualize_side_chains( core::pose::Pose & pose );

utility::vector1< core::Size >
get_all_residues( core::pose::Pose const & pose );

core::Size
get_unique_connection_res( core::pose::Pose const & pose, utility::vector1< core::Size > const & moving_partition_res );

} //modeler
} //stepwise
} //protocols

#endif
