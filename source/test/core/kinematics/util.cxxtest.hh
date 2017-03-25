// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/util.cxxtest.hh
/// @brief
/// @author Rocco Moretti (rmorettiase@gmail.com)
/// @author Steven Lewis smlewi@gmail.com test_utility_function_remodel_fold_tree_to_account_for_insertion

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/core/kinematics/utilities.hh>

// Unit Headers
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>

// Package Header
#include <core/kinematics/util.hh>

//Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/stream_util.hh>

// C++ Headers
#include <sstream>

using core::kinematics::FoldTree;
using core::kinematics::Edge;

static THREAD_LOCAL basic::Tracer TR( "core.kinematics.util.cxxtest.hh" );

class CoreKinematicsUtilTests : public CxxTest::TestSuite {

private:

	FoldTree mixedup_ft_;

public:

	void setUp()
	{
		core_init();

		//  +--------------------------2-----------------------------+
		//  |                                  +-1-+ +------3------+ |
		//  ^                                  v   ^ v             ^ v
		//  1 -> 2 -> 3 -> 4    5 <- 6 <- 7 <- 8    9    10    11   12 -> 13
		//  v                                            ^ v    ^
		//  +--------------------------5-----------------+ +--4-+
		//
		std::istringstream mixedup_ft_is("FOLD_TREE  EDGE 1 5 -1  EDGE 1 12 2 EDGE 12 13 -1 EDGE 12 9 3 EDGE 9 8 1 EDGE 8 6 -1 EDGE 1 10 5 EDGE 10 11 4");

		mixedup_ft_is >> mixedup_ft_;

		TR << "Mixed up foldtree: " << mixedup_ft_ << std::endl;
		runtime_assert( mixedup_ft_.connected() );
		runtime_assert( mixedup_ft_.check_fold_tree() );
	}

	void tearDown()
	{}

public:

	void test_jump_which_partitions() {

		using namespace core::kinematics;

		// --------------------------------- 1  2  3  4  5  6  7  8  9 10 11 12 13
		utility::vector1< bool > partition1{ 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1};
		utility::vector1< bool > partition2{ 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1};
		utility::vector1< bool > partition3{ 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0};
		utility::vector1< bool > partition4{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1};
		utility::vector1< bool > partition5{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0};

		utility::vector1< bool > badpartit1{ 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1}; // Different in peptide
		utility::vector1< bool > badpartit2{ 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0}; // Different in branches
		utility::vector1< bool > badpartit3{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0};
		utility::vector1< bool > badpartit4{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}; // all the same

		TS_ASSERT_EQUALS( jump_which_partitions( mixedup_ft_, partition1 ), 1 );
		TS_ASSERT_EQUALS( jump_which_partitions( mixedup_ft_, partition2 ), 2 );
		//TR << "Desired Partition 2: " << partition2 << std::endl;
		//TR << "Actual Partition 2: " << mixedup_ft_.partition_by_jump(2) << std::endl;
		TS_ASSERT_EQUALS( jump_which_partitions( mixedup_ft_, partition3 ), 3 );
		//TR << "Desired Partition 3: " << partition3 << std::endl;
		//TR << "Actual Partition 3: " << mixedup_ft_.partition_by_jump(3) << std::endl;
		TS_ASSERT_EQUALS( jump_which_partitions( mixedup_ft_, partition4 ), 4 );
		TS_ASSERT_EQUALS( jump_which_partitions( mixedup_ft_, partition5 ), 5 );
		TS_ASSERT_EQUALS( jump_which_partitions( mixedup_ft_, badpartit1 ), 0 );
		TS_ASSERT_EQUALS( jump_which_partitions( mixedup_ft_, badpartit2 ), 0 );
		TS_ASSERT_EQUALS( jump_which_partitions( mixedup_ft_, badpartit3 ), 0 );
		TS_ASSERT_EQUALS( jump_which_partitions( mixedup_ft_, badpartit4 ), 0 );

	}

	void test_get_foldtree_which_partitions() {

		using namespace core::kinematics;

		// --------------------------------- 1  2  3  4  5  6  7  8  9 10 11 12 13
		utility::vector1< bool > partition1{ 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1};
		utility::vector1< bool > partition2{ 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1};
		utility::vector1< bool > partition3{ 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0};
		utility::vector1< bool > partition4{ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1};
		utility::vector1< bool > partition5{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0};

		// If we pass in a partition which already corresponds with a jump, we should get back an identical FoldTree.
		TS_ASSERT( mixedup_ft_.is_equivalent( get_foldtree_which_partitions( mixedup_ft_, partition1 ) ) );
		TS_ASSERT( mixedup_ft_.is_equivalent( get_foldtree_which_partitions( mixedup_ft_, partition2 ) ) );
		TS_ASSERT( mixedup_ft_.is_equivalent( get_foldtree_which_partitions( mixedup_ft_, partition3 ) ) );
		TS_ASSERT( mixedup_ft_.is_equivalent( get_foldtree_which_partitions( mixedup_ft_, partition4 ) ) );
		TS_ASSERT( mixedup_ft_.is_equivalent( get_foldtree_which_partitions( mixedup_ft_, partition5 ) ) );

		// --------------------------------- 1  2  3  4  5  6  7  8  9 10 11 12 13
		utility::vector1< bool > newpartit1{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1}; // graft 6-8 to root
		utility::vector1< bool > newpartit2{ 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0}; // Split up some polymers
		utility::vector1< bool > newpartit3{ 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1}; // split only in polymer

		FoldTree new_fold_tree1( get_foldtree_which_partitions( mixedup_ft_, newpartit1 ) );
		TR << "FoldTree 1: " << new_fold_tree1 << std::endl;
		TR << "Desired Partition 1: " << newpartit1 << std::endl;
		TR << "Actual Partition 1:  " << new_fold_tree1.partition_by_jump(2) << std::endl;
		TS_ASSERT_FOLD_TREE_HAS_JUMP( new_fold_tree1, 1, 12 );
		TS_ASSERT_FOLD_TREE_HAS_JUMP( new_fold_tree1, 12, 9 );
		TS_ASSERT_FOLD_TREE_HAS_EDGE( new_fold_tree1, 12, 13 ); // Non-jump edge
		TS_ASSERT_FOLD_TREE_HAS_EDGE( new_fold_tree1, 8, 6 ); // Non-jump edge
		TS_ASSERT_FOLD_TREE_HAS_JUMP( new_fold_tree1, 1, 8 ); // A bit of white-box testing - feel free to change if behavior changes
		TS_ASSERT_EQUALS( jump_which_partitions( new_fold_tree1, newpartit1 ), 2 ); // white box

		FoldTree new_fold_tree2( get_foldtree_which_partitions( mixedup_ft_, newpartit2 ) );
		TR << "FoldTree 2: " << new_fold_tree2 << std::endl;
		TR << "Desired Partition 2: " << newpartit2 << std::endl;
		TR << "Actual Partition 2:  " << new_fold_tree2.partition_by_jump(5) << std::endl;
		TS_ASSERT_FOLD_TREE_HAS_JUMP( new_fold_tree2, 1, 10 );
		TS_ASSERT_FOLD_TREE_HAS_JUMP( new_fold_tree2, 10,11 );
		TS_ASSERT_FOLD_TREE_HAS_JUMP( new_fold_tree2, 2,  5 );
		TS_ASSERT_FOLD_TREE_HAS_EDGE( new_fold_tree2, 1, 2 ); //
		TS_ASSERT_FOLD_TREE_HAS_EDGE( new_fold_tree2, 3, 4 ); // Non-jump edge
		TS_ASSERT_FOLD_TREE_HAS_JUMP( new_fold_tree2, 10, 13 ); // A bit of white-box testing - feel free to change if behavior changes
		TS_ASSERT_FOLD_TREE_HAS_JUMP( new_fold_tree2, 10, 3 ); // white-box
		TS_ASSERT_EQUALS( jump_which_partitions( new_fold_tree2, newpartit2 ), 5 ); // white box

		FoldTree new_fold_tree3( get_foldtree_which_partitions( mixedup_ft_, newpartit3 ) );
		TR << "FoldTree 3: " << new_fold_tree3 << std::endl;
		TR << "Desired Partition 3: " << newpartit3 << std::endl;
		TR << "Actual Partition 3:  " << new_fold_tree3.partition_by_jump(6) << std::endl;
		TS_ASSERT_FOLD_TREE_HAS_JUMP( new_fold_tree3, 9, 8 );
		TS_ASSERT_FOLD_TREE_HAS_EDGE( new_fold_tree3, 7, 6 ); // Non-jump edge
		TS_ASSERT_FOLD_TREE_HAS_JUMP( new_fold_tree3, 1, 7 ); // A bit of white-box testing - feel free to change if behavior changes
		TS_ASSERT_EQUALS( jump_which_partitions( new_fold_tree3, newpartit3 ), 6 ); // white box, new jump for this purpose

	}

private:
};
