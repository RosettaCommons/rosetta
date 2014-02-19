// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  Test the algorithms that drive the KinematicMover.
/// @author Kale Kundert

#ifndef INCLUDED_protocols_kinematic_closure_KicFoldTree_CXXTEST_HH
#define INCLUDED_protocols_kinematic_closure_KicFoldTree_CXXTEST_HH

#define private public
#define protected public

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/core/kinematics/utilities.hh>

// Unit headers
#include <protocols/kinematic_closure/utilities.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>

using namespace std;
using namespace core;
using core::pose::Pose;
using core::import_pose::pose_from_pdb;
using core::kinematics::FoldTree;
using protocols::loops::Loop;
using protocols::loops::Loops;

class KicFoldTreeTests : public CxxTest::TestSuite {

public:

	void setUp() { core_init(); }

	FoldTree make_test_case(string path, Size start, Size stop, Size cut=0) {
		Pose pose; pose_from_pdb(pose, path);
		Loop loop(start, stop, cut);
		protocols::kinematic_closure::setup_fold_tree(pose, loop);
		return pose.fold_tree();
	}

	void test_ordinary_tree() {
		FoldTree tree = make_test_case(
				"protocols/kinematic_closure/inputs/3cla.pdb", 170, 181, 181);

		TS_ASSERT_FOLD_TREE_HAS_EDGE(tree,   1, 169);
		TS_ASSERT_FOLD_TREE_HAS_EDGE(tree, 169, 182);
		TS_ASSERT_FOLD_TREE_HAS_EDGE(tree, 183, 213);
		TS_ASSERT_FOLD_TREE_HAS_JUMP(tree, 169, 183);
	}

	void test_tree_with_ligand() {
		FoldTree tree = make_test_case(
				"protocols/kinematic_closure/inputs/1exm.pdb", 289, 300, 300);

		TS_ASSERT_FOLD_TREE_HAS_EDGE(tree,   1, 288);
		TS_ASSERT_FOLD_TREE_HAS_EDGE(tree, 288, 301);
		TS_ASSERT_FOLD_TREE_HAS_EDGE(tree, 302, 403);
		TS_ASSERT_FOLD_TREE_HAS_JUMP(tree, 288, 302);
		TS_ASSERT_FOLD_TREE_HAS_JUMP(tree,   1, 404);
	}

	void test_tree_with_symmetry() {}
};

#endif



