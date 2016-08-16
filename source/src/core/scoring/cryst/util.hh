// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_cryst_util_hh
#define INCLUDED_core_scoring_cryst_util_hh

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <iostream>


namespace core {
namespace scoring {
namespace cryst {

/// @brief fix the hydrogen bfactors in the pose
void fix_bfactorsH( core::pose::Pose & );
void fix_bfactorsMissing( core::pose::Pose & pose );


} // namespace constraints
} // namespace scoring
} // namespace core

#endif
