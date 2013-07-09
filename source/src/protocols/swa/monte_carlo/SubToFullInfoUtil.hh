// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SubToFullInfo.hh
/// @brief  Helper functions for SubToFull object, cached inside pose.
/// @author Rhiju Das

#ifndef INCLUDED_protocols_swa_monte_carlo_SubToFullInfoUtil_hh
#define INCLUDED_protocols_swa_monte_carlo_SubToFullInfoUtil_hh

#include <core/types.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/swa/monte_carlo/SubToFullInfo.fwd.hh>

namespace protocols {
namespace swa {
namespace monte_carlo {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of RNA centroid information inside the pose.

// Undefined, commenting out to fix PyRosetta build  SubToFullInfo const & sub_to_full_info_from_pose( core::pose::Pose const & pose );

SubToFullInfo &
nonconst_sub_to_full_info_from_pose( core::pose::Pose & pose );

void
reorder_sub_to_full_info_after_delete( core::pose::Pose & pose, core::Size const res_to_delete );

void
reorder_sub_to_full_info_after_append( core::pose::Pose & pose, core::Size const res_to_add );

void
reorder_sub_to_full_info_after_prepend( core::pose::Pose & pose, core::Size const res_to_add );

void
update_pdb_info_from_sub_to_full( core::pose::Pose & pose );


}
}
}
#endif
