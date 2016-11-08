// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/chainbreak_util.hh
/// @brief  Utility functions for scoring chainbreaks.
/// @author James Thompson

#ifndef INCLUDED_core_scoring_methods_chainbreak_util_hh
#define INCLUDED_core_scoring_methods_chainbreak_util_hh

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>

#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/conformation/Residue.fwd.hh>


namespace core {
namespace scoring {
namespace methods {

bool is_lower_cutpoint(
	core::Size residue,
	core::pose::Pose const & pose
);

bool is_upper_cutpoint(
	core::Size residue,
	core::pose::Pose const & pose
);

void find_cutpoint_variants(
	core::pose::Pose const & pose,
	core::kinematics::FoldTree const & tree,
	utility::vector1<int> & cutpoints
);

/// @brief helper function for looking at residue connections to get lower/upper partners
bool
lower_upper_connected_across_cutpoint( core::conformation::Residue const & lower_rsd,
																			 core::conformation::Residue const & upper_rsd );

/// @brief Instead of assuming cutpoint partner is simply cutpoint+1, find which residue connects via lower/upper.
core::Size
get_upper_cutpoint_partner_for_lower( core::pose::Pose const & pose, core::Size const lower_res );

/// @brief Instead of assuming cutpoint partner is simply res-1, find which residue connects via lower/upper.
core::Size
get_lower_cutpoint_partner_for_upper( core::pose::Pose const & pose, core::Size const upper_res );

} // namespace methods
} // namespace scoring
} // namespace core

#endif
