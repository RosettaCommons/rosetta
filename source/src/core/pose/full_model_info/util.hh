// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
	update_pdb_info_from_full_model_info( core::pose::Pose & pose );

	utility::vector1< char >
	figure_out_conventional_chains_from_full_model_info( pose::Pose const & pose );

	utility::vector1< Size >
	figure_out_chain_numbers_from_full_model_info( pose::Pose & pose );

	utility::vector1< Size >
	figure_out_chain_numbers_from_full_model_info_const( pose::Pose const & pose );

	utility::vector1< Size >
	get_chains_full( pose::Pose const & pose );

	bool
	check_full_model_info_OK( pose::Pose const & pose );

	utility::vector1< Size > const &
	get_res_list_from_full_model_info( pose::Pose & pose );

	utility::vector1< Size > const &
	get_res_list_from_full_model_info_const( pose::Pose const & pose );

	utility::vector1< utility::vector1< Size > >
	get_move_elements_from_full_model_info( pose::Pose & pose );

	utility::vector1< utility::vector1< Size > >
	get_move_elements_from_full_model_info_const( pose::Pose const & pose );

	utility::vector1< Size >
	get_moving_res_from_full_model_info( pose::Pose & pose );

	utility::vector1< Size >
	get_fixed_domain_from_full_model_info_const( pose::Pose const & pose );

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

	// Undefined, commenting out to fix PyRosetta build
	// utility::vector1< Size >
	// figure_out_pose_domain_map_const( core::pose::Pose const & pose );

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

	Size
	get_number_missing_residue_connections( pose::Pose & pose );


} //full_model_info
} //pose
} //core
#endif
