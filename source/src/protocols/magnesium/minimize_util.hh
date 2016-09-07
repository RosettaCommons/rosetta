// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/magnesium/minimize_util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_magnesium_minimize_util_HH
#define INCLUDED_protocols_magnesium_minimize_util_HH

#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

namespace protocols {
namespace magnesium {

void
minimize_magnesium_and_hydration_shell( core::pose::Pose & pose,
	utility::vector1< core::Size > const & mg_res,
	core::scoring::ScoreFunctionCOP minimize_scorefxn = nullptr,
	core::Distance const mg_coord_cst_dist = 0.2 );

void
minimize_magnesium_and_hydration_shell( core::pose::Pose & pose /*for viewing*/,
	utility::vector1< core::pose::PoseOP > & pose_list,
	utility::vector1< core::Size > const & mg_res,
	core::scoring::ScoreFunctionCOP minimize_scorefxn = nullptr,
	core::Distance const mg_coord_cst_dist = 0.2 );

void
update_mg_hoh_fold_tree( core::pose::Pose & pose );

core::id::AtomID
get_closest_non_hoh_contact( core::pose::Pose const & pose, core::Size const i, std::string const & exclude_rsd = "" );

} //magnesium
} //protocols



#endif
