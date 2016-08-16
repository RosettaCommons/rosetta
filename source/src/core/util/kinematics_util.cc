// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/util/kinematics_util.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <core/util/kinematics_util.hh>

// Project headers
#include <core/chemical/VariantType.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

//Auto Headers
#include <utility/vector1.hh>

namespace core {
namespace util {

void remove_cutpoint_variants(core::pose::Pose & pose) {
	const core::kinematics::FoldTree& tree(pose.fold_tree());
	for ( unsigned i = 1; i <= pose.total_residue(); ++i ) {
		if ( !tree.is_cutpoint(i) || i >= (pose.total_residue() - 1) ) {
			continue;
		}

		core::pose::remove_variant_type_from_pose_residue(pose, core::chemical::CUTPOINT_LOWER, i);
		core::pose::remove_variant_type_from_pose_residue(pose, core::chemical::CUTPOINT_UPPER, i + 1);
	}
}

}  // namespace core
}  // namespace util
