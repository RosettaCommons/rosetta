// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/ConsecutiveTreeBuilder.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/ConsecutiveTreeBuilder.hh>

// Package headers
#include <protocols/nonlocal/CutFinder.hh>
#include <protocols/nonlocal/NLFragmentGroup.hh>
#include <protocols/nonlocal/NLGrouping.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace nonlocal {

void ConsecutiveTreeBuilder::build(const NLGrouping& grouping,
                                   core::pose::Pose* pose,
                                   core::fragment::SecondaryStructureOP ss) {
  using core::Size;
  using core::kinematics::FoldTree;

  FoldTree tree(pose->fold_tree());
  for (Size i = 1; i <= grouping.num_groups() - 1; ++i) {
    const NLFragmentGroup& fragment_1 = grouping.groups(i);
    const NLFragmentGroup& fragment_2 = grouping.groups(i+1);
    Size cutpoint = CutFinder::choose_cutpoint(fragment_1.stop() + 1, fragment_2.start() - 1, ss);

    // Identify the central residue of each fragment
    Size midpoint_1 = (fragment_1.start() + fragment_1.stop()) / 2;
    Size midpoint_2 = (fragment_2.start() + fragment_2.stop()) / 2;
    int jump_id = tree.new_jump(midpoint_1, midpoint_2, cutpoint);
    tree.set_jump_atoms(jump_id, "CA", "CA");
  }

  // Update the pose's FoldTree
  pose->fold_tree(tree);
}

}  // namespace nonlocal
}  // namespace protocols
