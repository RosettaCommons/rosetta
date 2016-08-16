// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/ClusteringTreeNodeTest.cxxtest.hh
/// @brief
/// @author Domini Gront

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <numeric/ClusteringTreeNode.hh>

#include <core/types.hh>
#include <iostream>

using namespace numeric;

class ClusteringTreeNodeTest : public CxxTest::TestSuite {
public:

	ClusteringTreeNodeTest() {};

	// Shared initialization goes here.
	void setUp() {
	}

	void test_node_walk() {

	    core::Size result[] = {1,2,4,5,3};
	    core::Size result_leaves[] = {4,5,3};
	    ClusteringTreeNodeOP lower_left_left = ClusteringTreeNode::newClusteringTreeNode(4);
	    ClusteringTreeNodeOP lower_left_right = ClusteringTreeNode::newClusteringTreeNode(5);
	    ClusteringTreeNodeOP left = ClusteringTreeNode::newClusteringTreeNode(2,lower_left_left,lower_left_right);
	    ClusteringTreeNodeOP right = ClusteringTreeNode::newClusteringTreeNode(3);
	    ClusteringTreeNodeOP root = ClusteringTreeNode::newClusteringTreeNode(1,left,right);

    	    utility::vector1<Size> ids;
	    root->copy_member_ids(ids);
	    assert( ids.size() == 5 );
	    for(Size i=1;i<=ids.size();i++) {
		TS_ASSERT( ids[i] == result[i-1] );
	    }
	    //std::cout <<std::endl;
	    ids.clear();
	    root->copy_leaf_ids(ids);
	    assert( ids.size() == 3 );
	    for(Size i=1;i<=ids.size();i++) {
		TS_ASSERT( ids[i] == result_leaves[i-1] );
	    }

	}

	// Shared finalization goes here.
	void tearDown() {}
};
