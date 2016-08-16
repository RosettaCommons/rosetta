// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/util/disulfide_util.hh
/// @brief A collection of procedures for manipulating disulfide bonds
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date 4/30/2009

#ifndef INCLUDED_core_util_disulfide_util_HH
#define INCLUDED_core_util_disulfide_util_HH

// Project Headers
#include <core/types.hh>

#include <core/kinematics/MoveMap.fwd.hh>

#include <core/pack/task/PackerTask.fwd.hh>

#include <core/pose/Pose.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/vector1.fwd.hh>

// C++ headers

#include <utility/vector1.hh>


namespace core {
namespace util {

/// @brief Rebuild a pair of cysteines (and possibly surrounding residues) so
///  that they form a near-ideal disulfide bond.  Supports symmetric poses.
void rebuild_disulfide( core::pose::Pose & pose,
	core::Size lower_res, core::Size upper_res,
	core::pack::task::PackerTaskOP packer_task = 0,
	core::scoring::ScoreFunctionOP packer_score = 0,
	core::kinematics::MoveMapOP mm = 0,
	core::scoring::ScoreFunctionOP minimizer_score = 0 );

/// @brief Rebuild a number of pairs of cysteines (and possibly surrounding
///  residues) so that they form near-ideal disulfide bonds
void rebuild_disulfide( core::pose::Pose & pose,
	utility::vector1<std::pair<core::Size, core::Size> > disulfides,
	core::pack::task::PackerTaskOP packer_task = 0,
	core::scoring::ScoreFunctionOP packer_score = 0,
	core::kinematics::MoveMapOP mm = 0,
	core::scoring::ScoreFunctionOP minimizer_score = 0 );

} // util
} // core

#endif //INCLUDED_protocols_toolbox_disulfide_util_HH

