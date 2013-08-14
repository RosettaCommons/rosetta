// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/LoopsUtil.cxxtest.hh
/// @brief test suite for protocols/loops/util
/// @author Christopher Miles (cmiles@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Package headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>

// Project headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>

#include <iostream>

namespace {

using core::kinematics::Edge;
using core::kinematics::FoldTree;
using core::pose::Pose;
using core::pose::PoseOP;
using protocols::loops::Loop;
using protocols::loops::Loops;
using protocols::loops::fold_tree_from_loops;

class LoopsFoldTreeTest : public CxxTest::TestSuite {

private:
	PoseOP pose_;
	PoseOP pose_multichain_;

public:
  void setUp() {
    protocols_init();
		pose_ = core::import_pose::pose_from_pdb("protocols/loops/2GB3.pdb");
		pose_multichain_ = core::import_pose::pose_from_pdb("protocols/loops/4DZM.pdb");
  }

	void test_SingleLoopFoldTree() {
		FoldTree ft(pose_->total_residue());

		Loops loops;
		loops.add_loop(Loop(5, 12, 8));

		// generate a FoldTree from the Loops instance
		fold_tree_from_loops(*pose_, loops, ft);

		// We should have the following FoldTree:
		// 1 --- 4 --- 8 9 -- 13 -------- N
		// ----->*----->|<-----*---------->
		//       |             |
		//       |_____________|
		//

		// This test will check for each part of the FoldTree in turn:
		// 1) Edge from  1 -> 4
		TS_ASSERT(ft.edge_label(1, 4) == Edge::PEPTIDE);

		// 2) Edge from  4 -> 8
		TS_ASSERT(ft.edge_label(4, 8) == Edge::PEPTIDE);

		// 3) Edge from 13 -> 9
		TS_ASSERT(ft.edge_label(13, 9) == Edge::PEPTIDE);

		// 4) Edge from 13 -> N
		TS_ASSERT(ft.edge_label(13, pose_->total_residue()) == Edge::PEPTIDE);

		// 5) Jump from  4 -> 13
		TS_ASSERT(ft.edge_label(4, 13) == 1);
	}

	void test_MultiLoopFoldTree() {
		FoldTree ft(pose_->total_residue());

		Loops loops;
		loops.add_loop(Loop(5, 12, 8));
		loops.add_loop(Loop(20, 28, 24));

		// generate a FoldTree from the Loops instance
		fold_tree_from_loops(*pose_, loops, ft);

		// We should have the following FoldTree:
		// 1 --- 4 --- 8 9 -- 13 ---- 19 --- 24 25 --- 29 ----- N
		// ----->*----->|<-----*------>*------>|<------*-------->
		//       |             |       |               |
		//       |_____________|       |_______________|
		//

		// Test first loop
		// 1) Edge from  1 -> 4
		TS_ASSERT(ft.edge_label(1, 4) == Edge::PEPTIDE);

		// 2) Edge from  4 -> 8
		TS_ASSERT(ft.edge_label(4, 8) == Edge::PEPTIDE);

		// 3) Edge from 13 -> 9
		TS_ASSERT(ft.edge_label(13, 9) == Edge::PEPTIDE);

		// 4) Edge from 13 -> 19
		TS_ASSERT(ft.edge_label(13, 19) == Edge::PEPTIDE);

		// 5) Jump from  4 -> 13
		TS_ASSERT(ft.edge_label(4, 13) == 1);

		// Test second loop
		// 1) Edge from  19 -> 24
		TS_ASSERT(ft.edge_label(19, 24) == Edge::PEPTIDE);

		// 2) Edge from  29 -> 25
		TS_ASSERT(ft.edge_label(29, 25) == Edge::PEPTIDE);

		// 3) Edge from 29 -> N
		TS_ASSERT(ft.edge_label(29, pose_->total_residue()) == Edge::PEPTIDE);

		// 4) Jump from 19 -> 29
		TS_ASSERT(ft.edge_label(19, 29) == 2);
	}

	void test_NTerminalLoopFoldTree() {
		FoldTree ft(pose_->total_residue());

		Loops loops;
		loops.add_loop(Loop(1, 12, 8));

		// generate a FoldTree from the Loops instance
		fold_tree_from_loops(*pose_, loops, ft);

		// We should have the following FoldTree:
		// 1 ---------------- 13 ------------------------------ N
		// <------------------*--------------------------------->
		//

		// Test N-term loop
		// 1) Edge from  13 -> 1
		TS_ASSERT(ft.edge_label(13, 1) == Edge::PEPTIDE);

		// 2) Edge from  13 -> N
		TS_ASSERT(ft.edge_label(13, pose_->total_residue()) == Edge::PEPTIDE);

	}

	void test_NTerminalLoopFoldTreeWithStupidCutpoint() {
		FoldTree ft(pose_->total_residue());

		Loops loops;
		loops.add_loop(Loop(1, 12, 8));

		// generate a FoldTree from the Loops instance
		fold_tree_from_loops(*pose_, loops, ft, true);

		// We should have the following FoldTree:
		// 1 --------- 8 9 -- 13 -------- N
		// *----------->|<-----*---------->
		// |                   |
		// |___________________|
		//

		// Test N-term loop
		// 1) Edge from  1 -> 8
		TS_ASSERT(ft.edge_label(1, 8) == Edge::PEPTIDE);

		// 2) Edge from  13 -> 9
		TS_ASSERT(ft.edge_label(13, 9) == Edge::PEPTIDE);

		// 3) Edge from  13 -> N
		TS_ASSERT(ft.edge_label(13, pose_->total_residue()) == Edge::PEPTIDE);

		// 4) Jump from 1 -> 13
		TS_ASSERT(ft.edge_label(1, 13) == 1);
	}

	void test_CTerminalLoopFoldTree() {
		FoldTree ft(pose_->total_residue());

		Loops loops;
		loops.add_loop(Loop(30, pose_->total_residue(), 35));

		// generate a FoldTree from the Loops instance
		fold_tree_from_loops(*pose_, loops, ft);

		// We should have the following FoldTree:
		// 1 ---------------- 29 ------------------------------ N
		// -------------------*--------------------------------->
		//

		// Test C-term loop
		// 1) Edge from  1 -> 29
		TS_ASSERT(ft.edge_label(1, 29) == Edge::PEPTIDE);

		// 2) Edge from  29 -> N
		TS_ASSERT(ft.edge_label(29, pose_->total_residue()) == Edge::PEPTIDE);
	}

	void test_CTerminalLoopFoldTreeWithStupidCutpoint() {
		FoldTree ft(pose_->total_residue());

		Loops loops;
		loops.add_loop(Loop(30, pose_->total_residue(), 35));

		// generate a FoldTree from the Loops instance
		fold_tree_from_loops(*pose_, loops, ft, true);

		// We should have the following FoldTree:
		// 1 ---------- 29 --- 35 36 -------- N
		// ------------>*------->|<-----------*
		//              |                     |
		//              |_____________________|
		//

		// Test N-term loop
		// 1) Edge from  1 -> 29
		TS_ASSERT(ft.edge_label(1, 29) == Edge::PEPTIDE);

		// 2) Edge from  29 -> 35
		TS_ASSERT(ft.edge_label(29, 35) == Edge::PEPTIDE);

		// 3) Edge from  N -> 36
		TS_ASSERT(ft.edge_label(pose_->total_residue(), 36) == Edge::PEPTIDE);

		// 4) Jump from 29 -> N
		TS_ASSERT(ft.edge_label(29, pose_->total_residue()) == 1);
	}

	void test_Multichain_CtermLoop() {
		FoldTree ft(pose_multichain_->total_residue());

		// pose_multichain_ has two 31-residue chains

		Loops loops;
		loops.add_loop(Loop(25, 31, 31));

		// generate a FoldTree from the Loops instance
		fold_tree_from_loops(*pose_multichain_, loops, ft);

		// Test N-term loop
		// 1) Edge from  1 -> 24
		TS_ASSERT(ft.edge_label(1, 24) == Edge::PEPTIDE);

		// 2) Edge from  24 -> 31
		TS_ASSERT(ft.edge_label(24, 31) == Edge::PEPTIDE);

		// 3) Edge from  32 -> 62
		TS_ASSERT(ft.edge_label(32, 62) == Edge::PEPTIDE);

		// 4) Jump from 24 -> 32
		TS_ASSERT(ft.edge_label(24, 32) == 1);
	}

	void test_Multichain_NtermLoop() {
		FoldTree ft(pose_multichain_->total_residue());

		// pose_multichain_ has two 31-residue chains
		Loops loops;
		loops.add_loop(Loop(32, 38, 32));

		// generate a FoldTree from the Loops instance
		fold_tree_from_loops(*pose_multichain_, loops, ft);

		TS_ASSERT(ft.edge_label(1, 31) == Edge::PEPTIDE);

		TS_ASSERT(ft.edge_label(39, 32) == Edge::PEPTIDE);

		TS_ASSERT(ft.edge_label(39, 62) == Edge::PEPTIDE);

		TS_ASSERT(ft.edge_label(31, 39) == 1);
	}

	void test_Multichain_NandCtermLoop() {
		FoldTree ft(pose_multichain_->total_residue());

		// pose_multichain_ has two 31-residue chains
		Loops loops;
		loops.add_loop(Loop(25, 31, 31));
		loops.add_loop(Loop(32, 38, 32));

		// generate a FoldTree from the Loops instance
		fold_tree_from_loops(*pose_multichain_, loops, ft);

		TS_ASSERT(ft.edge_label(1, 24) == Edge::PEPTIDE);

		TS_ASSERT(ft.edge_label(24, 31) == Edge::PEPTIDE);

		TS_ASSERT(ft.edge_label(39, 32) == Edge::PEPTIDE);

		TS_ASSERT(ft.edge_label(39, 62) == Edge::PEPTIDE);

		TS_ASSERT(ft.edge_label(24, 39) == 1);

	}
};

}  // anonymous namespace
