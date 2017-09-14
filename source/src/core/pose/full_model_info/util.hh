// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/FullModelInfo.hh
/// @brief  Helper functions for SubToFull object, cached inside pose.
/// @author Rhiju Das

#ifndef INCLUDED_core_pose_full_model_info_FullModelInfoUtil_hh
#define INCLUDED_core_pose_full_model_info_FullModelInfoUtil_hh

#include <core/types.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/full_model_info/FullModelInfo.fwd.hh>
#include <core/pose/full_model_info/FullModelParameters.fwd.hh>
#include <map>

namespace core {
namespace pose {
namespace full_model_info {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of RNA centroid information inside the pose.

void
reorder_full_model_info_after_delete( core::pose::Pose & pose, core::Size const res_to_delete );

void
reorder_full_model_info_after_append( core::pose::Pose & pose, core::Size const res_to_add, Size const offset = 1 );

void
reorder_full_model_info_after_prepend( core::pose::Pose & pose, core::Size const res_to_add, Size const offset = 1 );

void
update_res_list_in_full_model_info_and_pdb_info( pose::Pose & pose, utility::vector1< Size > const & res_list_new );

void
update_pose_objects_from_full_model_info( core::pose::Pose & pose );

void
update_pdb_info_from_full_model_info( core::pose::Pose & pose );

void
update_constraint_set_from_full_model_info( core::pose::Pose & pose );

void
update_disulfides_from_full_model_info( pose::Pose & pose );

utility::vector1< char >
figure_out_conventional_chains_from_full_model_info( pose::Pose const & pose );

utility::vector1< Size >
figure_out_chain_numbers_from_full_model_info( pose::Pose & pose );

utility::vector1< Size >
figure_out_chain_numbers_from_full_model_info_const( pose::Pose const & pose );

utility::vector1< Size >
figure_out_dock_domain_map_from_full_model_info_const( pose::Pose const & pose );

utility::vector1< Size >
get_chains_full( pose::Pose const & pose );

utility::vector1< Size >
get_chains_from_cutpoint_open( utility::vector1< Size > const & cutpoint_open, Size const nres );

utility::vector1< Size >
get_cutpoint_open_from_chains( utility::vector1< Size > const & chains );

utility::vector1< Size >
get_sample_res_for_pose( pose::Pose const & pose );

bool
check_full_model_info_OK( pose::Pose const & pose );

utility::vector1< Size > const &
get_res_list_from_full_model_info( pose::Pose & pose );

utility::vector1< Size > const &
get_res_list_from_full_model_info_const( pose::Pose const & pose );

utility::vector1< Size >
get_res_list_const( pose::Pose const & pose );

utility::vector1< utility::vector1< Size > >
get_move_elements_from_full_model_info( pose::Pose & pose );

utility::vector1< utility::vector1< Size > >
get_move_elements_from_full_model_info_const( pose::Pose const & pose );

utility::vector1< Size >
get_moving_res_from_full_model_info( pose::Pose & pose );

utility::vector1< Size >
get_moving_res_from_full_model_info_const( pose::Pose const & pose );

utility::vector1< Size >
get_fixed_domain_from_full_model_info_const( pose::Pose const & pose );

utility::vector1< Size >
get_input_domain_from_full_model_info_const( pose::Pose const & pose );

core::Size
sub_to_full( core::Size const i, core::pose::Pose const & pose );

utility::vector1< Size >
sub_to_full( utility::vector1< Size > const & res_list, core::pose::Pose const & pose );

core::Size
full_to_sub( core::Size const i, core::pose::Pose const & pose );

utility::vector1< Size >
full_to_sub( utility::vector1< Size > const & res_list, core::pose::Pose const & pose );

Size
full_model_size( pose::Pose & pose );

utility::vector1< Size >
figure_out_pose_domain_map( core::pose::Pose & pose );

utility::vector1< Size >
figure_out_pose_domain_map_const( core::pose::Pose const & pose );

core::conformation::Residue const &
get_residue( Size const seqpos_in_full_model,
	pose::Pose const & pose );

utility::vector1< int >
get_res_num_from_pdb_info( pose::Pose const & pose );

utility::vector1< char >
get_chains_from_pdb_info( pose::Pose const & pose );

Size
get_chain_for_full_model_resnum( Size const & resnum, pose::Pose const & pose );

Size
get_chain_for_resnum( Size const & resnum, pose::Pose const & pose );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// deprecate soon in favor of direct loop-graph call
//
/// @brief Finds the number of missing residues in a pose.
/// @details This function returns the number of missing residues in the pose.
/// The pose is passed by nonconst reference, so that the full_model_info can be
/// setup, if needed.
Size
get_number_missing_residues_and_connections( pose::Pose & pose );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// deprecate soon in favor of direct loop-graph call
//
/// @brief Finds the number of missing residues in a pose.
/// @details This function returns the number of missing residues in the pose.
/// The missing_residues vector is passed by nonconst reference, so its values
/// can be modified and accessed by this function and the calling method.
Size
get_number_missing_residues_and_connections( pose::Pose const & pose,
	utility::vector1< char > & missing_residues );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// deprecate soon in favor of direct loop-graph call
//
/// @brief Finds the number of missing residues in a pose.
/// @details This function returns the number of missing residues in the pose.
/// The loop_suites vector is passed by nonconst reference, so its values
/// can be modified and accessed by this function and the calling method.
Size
get_number_missing_residues_and_connections( pose::Pose const & pose,
	utility::vector1< utility::vector1< Size > > & loop_suites );

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// sort of a hodge-podge function -- might be better to unify with LoopGraph?
//
/// @brief Finds the number of missing residues in a pose.
/// @details This function returns the number of missing residues in the pose.
/// The missing_res and loop_suites vectors are passed by nonconst reference, so their
/// values can be modified and accessed by this function and the calling method.
Size
get_number_missing_residues_and_connections( pose::Pose const & pose,
	utility::vector1< char > & missing_residues,
	utility::vector1< utility::vector1< Size > > & loop_suites );

bool
check_all_residues_sampled( pose::Pose const & pose );

std::map< std::pair< Size, Size >, std::pair< Size, Size > >
get_preferred_jump_pair_for_docking_domains( FullModelInfo const & full_model_info );

utility::vector1< std::pair< core::Size, core::Size > >
get_chain_connections( core::pose::Pose const & pose );

utility::vector1< core::Size >
get_connection_domains( utility::vector1< std::pair< core::Size, core::Size > > const & chain_connections,
	core::Size const nchains );

bool
check_sample_sugar_in_full_model_info( core::pose::Pose const & pose,
	core::Size const i );

void
append_virtual_residue_to_full_model_info( core::pose::Pose & pose );




} //full_model_info
} //pose
} //core
#endif
