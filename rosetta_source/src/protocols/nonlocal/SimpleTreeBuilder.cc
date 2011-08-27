// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/SimpleTreeBuilder.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/SimpleTreeBuilder.hh>

// Package headers
#include <protocols/nonlocal/NLGrouping.hh>

// Project headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace nonlocal {

void SimpleTreeBuilder::set_up(const NLGrouping& grouping, core::pose::Pose* pose) {
  core::kinematics::FoldTree simple_tree(pose->total_residue());
  pose->fold_tree(simple_tree);
}

}  // namespace nonlocal
}  // namespace protocols
