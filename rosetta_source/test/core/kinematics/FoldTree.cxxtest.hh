// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/FoldTree.cxxtest.hh
/// @brief  test suite for core::kinematics::FoldTree.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Test Headers
#include <cxxtest/TestSuite.h>

// Unit Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>

// Project headers
#include <test/core/init_util.hh>

using core::kinematics::FoldTree;
using core::pose::Pose;

class FoldTreeTest : public CxxTest::TestSuite {
 public:
  Pose pose_;

  void setUp() {
    core_init();
    core::import_pose::pose_from_pdb(pose_, "core/kinematics/test.pdb");
  }

  // Ensure that repeated calls to hash value yield identical results when
  // applied to the same FoldTree
  void test_hash_value_unmodified() {
    TS_ASSERT_EQUALS(pose_.fold_tree().hash_value(),
                     pose_.fold_tree().hash_value());
  }

  // Ensure that operations that modify the FoldTree have an effect on the
  // calculation of the hash value
  void test_hash_value_modified() {
    FoldTree mod_tree(pose_.fold_tree());
    mod_tree.new_jump(1, 3, 2);

    size_t hash_orig = pose_.fold_tree().hash_value();
    size_t hash_mod = mod_tree.hash_value();
    TS_ASSERT_DIFFERS(hash_orig, hash_mod);
  }

  /// @brief A quick check to see whether the FoldTree considers the final
  /// residues in a simple fold tree to be a cutpoint
  void test_last_residue_is_cutpoint() {
    const FoldTree& tree(pose_.fold_tree());
    TS_ASSERT(tree.is_cutpoint(pose_.total_residue()));
  }

  void test_boundary_right() {
    const FoldTree& tree = pose_.fold_tree();
    const unsigned n = pose_.total_residue();

    for (unsigned i = 1; i <= pose_.total_residue(); ++i) {
      if (tree.is_root(i))
        continue;

      // Latter is a shortcut to the former
      TS_ASSERT_EQUALS(n, tree.get_residue_edge(i).stop());
      TS_ASSERT_EQUALS(n, tree.boundary_right(i));
    }
  }

  void test_boundary_left() {
    const FoldTree& tree = pose_.fold_tree();

    for (unsigned i = 1; i < pose_.total_residue(); ++i) {
      if (tree.is_root(i) || i == tree.nres())
        continue;

      // Latter is a shortcut to the former
      TS_ASSERT_EQUALS(1, tree.get_residue_edge(i).start());
      TS_ASSERT_EQUALS(1, tree.boundary_left(i));
    }
  }
};
