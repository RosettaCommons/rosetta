// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_core_util_cryst_util_hh
#define INCLUDED_core_util_cryst_util_hh

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <iostream>


namespace core {
namespace util {

///@brief get the optimal weight on the xtal ML energy
core::Real getMLweight( core::scoring::ScoreFunction & scorefxn, core::pose::Pose & );

///@brief get the optimal weight on the xtal ML energy considering only movable DOFs
core::Real getMLweight( core::scoring::ScoreFunction & scorefxn,  core::pose::Pose & , core::kinematics::MoveMap &);

///@brief get the optimal weight on the xtal ML energy
core::Real getMLweight_cart( core::scoring::ScoreFunction & scorefxn, core::pose::Pose & );

///@brief get the optimal weight on the xtal ML energy considering only movable DOFs
core::Real getMLweight_cart( core::scoring::ScoreFunction & scorefxn, core::pose::Pose & , core::kinematics::MoveMap &);


}
} // namespace core

#endif
