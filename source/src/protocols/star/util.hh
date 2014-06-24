// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/star/util.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_STAR_UTIL_HH
#define INCLUDED_PROTOCOLS_STAR_UTIL_HH

// C/C++ headers
#include <string>

// Project headers
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace star {

/// @brief Writes pose to disk with the specified filename if -abinitio:debug is enabled
void emit_intermediate(const core::pose::Pose& pose, const std::string& silent_filename);

/// @brief Restores simple kinematics to pose
void simple_fold_tree(core::pose::Pose & pose);

/// @brief Converts pose to centroid residue type set
void to_centroid(core::pose::Pose & pose);

}  // namespace star
}  // namespace protocols

#endif  // PROTOCOLS_STAR_UTIL_HH_
