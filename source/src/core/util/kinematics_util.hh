// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/util/kinematics_util.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_CORE_UTIL_KINEMATICS_UTIL_HH
#define INCLUDED_CORE_UTIL_KINEMATICS_UTIL_HH

// Project headers
#include <core/pose/Pose.fwd.hh>

//Auto Headers
#include <utility/vector1.hh>
namespace core {
namespace util {

// @brief Returns true if `pose` has a chainbreak, false otherwise
// Undefinded, commenting out to fix PyRosetta Buuild  bool has_chainbreak(const core::pose::Pose& pose);

// @brief Adds cutpoint variants to `pose` by scanning the fold tree
// Use core::pose::correctly_add_cutpoint_variants() instead.
//void add_cutpoint_variants(core::pose::Pose & pose);

/// @brief Removes cutpoint variants from `pose` by scanning the fold tree
void remove_cutpoint_variants(core::pose::Pose & pose);

}  // namespace core
}  // namespace util

#endif  // CORE_UTIL_KINEMATICS_UTIL_HH_
