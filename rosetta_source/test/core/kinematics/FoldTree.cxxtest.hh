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
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Project headers
#include <test/core/init_util.hh>

// C/C++
#include <iostream>

using core::kinematics::FoldTree;
using core::pose::Pose;

namespace {

basic::Tracer TR("core.kinematics.FoldTree.cxxtest");

class FoldTreeTest : public CxxTest::TestSuite {
 public:
  Pose pose_;

  void setUp() {
    core_init();
    core::import_pose::pose_from_pdb(pose_, "core/kinematics/test.pdb");
  }

  void tearDown() {}

  // Ensure that repeated calls to hash value yield identical results when
	// applied to the same FoldTree
  void test_hash_value_unmodified() {
		TS_ASSERT_EQUALS(pose_.fold_tree().hash_value(),
										 pose_.fold_tree().hash_value());
  }

  // Ensure that operations that modify the FoldTree have an effect on the
  // calculation of the hash value
  void test_hash_value_modified() {
		using std::endl;

    FoldTree mod_tree(pose_.fold_tree());
    mod_tree.new_jump(1, 3, 2);

    size_t hash_orig = pose_.fold_tree().hash_value();
    size_t hash_mod = mod_tree.hash_value();

    TR << "Original fold-tree: " << hash_orig << endl
       << "Modified fold-tree: " << hash_mod << endl;

    TS_ASSERT_DIFFERS(hash_orig, hash_mod);
  }

	/// @brief A quick check to see whether the FoldTree considers the final
	/// residues in a simple fold tree to be a cutpoint
	void test_last_residue_is_cutpoint() {
		FoldTree tree(pose_.fold_tree());
		TS_ASSERT(tree.is_cutpoint(pose_.total_residue()));
	}
};
}  // anonymous namespace
