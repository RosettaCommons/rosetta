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
#include <core/pose/full_model_info/FullModelInfo.fwd.hh>

namespace core {
namespace pose {
namespace full_model_info {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of RNA centroid information inside the pose.

	void
	reorder_full_model_info_after_delete( core::pose::Pose & pose, core::Size const res_to_delete );

	void
	reorder_full_model_info_after_append( core::pose::Pose & pose, core::Size const res_to_add );

	void
	reorder_full_model_info_after_prepend( core::pose::Pose & pose, core::Size const res_to_add );

	void
	update_res_list_in_full_model_info_and_pdb_info( pose::Pose & pose, utility::vector1< Size > const & res_list_new );

	void
	update_pdb_info_from_full_model_info( core::pose::Pose & pose );

	utility::vector1< Size >
	figure_out_chains_from_full_model_info( pose::Pose & pose );

	void
	fill_full_model_info_from_command_line( pose::Pose & pose );

	void
	fill_full_model_info_from_command_line( utility::vector1< pose::PoseOP > pose_list );

	bool
	check_full_model_info_OK( pose::Pose const & pose );

	utility::vector1< Size > const &
	get_res_list_from_full_model_info( pose::Pose & pose );

	utility::vector1< Size > const &
	get_res_list_from_full_model_info_const( pose::Pose const & pose );

	utility::vector1< utility::vector1< Size > >
	get_move_elements_from_full_model_info( pose::Pose & pose );

	utility::vector1< Size >
	get_moving_res_from_full_model_info( pose::Pose & pose );

	core::Size
	sub_to_full( core::Size const i, core::pose::Pose & pose );

	core::Size
	full_to_sub( core::Size const i, core::pose::Pose & pose );

}
}
}
#endif
