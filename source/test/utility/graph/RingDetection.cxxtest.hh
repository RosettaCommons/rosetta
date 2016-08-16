// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/graph/ring_detection.cxxtest.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)

// This has to be first in order to set things up dor the ring_detection.hh header below.
#include <boost/graph/adjacency_list.hpp>

// Test headers
#include <cxxtest/TestSuite.h>
#include <utility/graph/ring_detection.hh>

#include <core/types.hh>

#include <utility/vector1.hh>

// Use this for simple integer VD of verticies.
typedef boost::adjacency_list< boost::listS, boost::vecS, boost::undirectedS > Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor VD;
typedef boost::graph_traits<Graph>::edge_descriptor ED;

namespace {

class RingDetectionTest : public CxxTest::TestSuite {
public:

	void setUp() {
	}

	void tearDown() {
	}

	void test_benzene() {
		Graph g(12);
		boost::add_edge(0, 1, g);
		boost::add_edge(1, 2, g);
		boost::add_edge(2, 3, g);
		boost::add_edge(3, 4, g);
		boost::add_edge(4, 5, g);
		boost::add_edge(5, 0, g);
		//hydrogens
		boost::add_edge(0, 6, g);
		boost::add_edge(1, 7, g);
		boost::add_edge(2, 8, g);
		boost::add_edge(3, 9, g);
		boost::add_edge(4, 10,g);
		boost::add_edge(5, 11,g);

		for( core::Size ii(0); ii < 6; ++ii ) {
			TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( ii, g ), 6);
		}
		for( core::Size ii(6); ii < 12; ++ii ) {
			TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( ii, g ), 999999);
		}
	} 

	void test_cubane() {
		Graph g(8);
		//face 1
		boost::add_edge(0, 1, g);
		boost::add_edge(1, 2, g);
		boost::add_edge(2, 3, g);
		boost::add_edge(3, 0, g);
		//face 2
		boost::add_edge(4, 5, g);
		boost::add_edge(5, 6, g);
		boost::add_edge(6, 7, g);
		boost::add_edge(7, 4, g);
		//connections
		boost::add_edge(0, 4, g);
		boost::add_edge(1, 5, g);
		boost::add_edge(2, 6, g);
		boost::add_edge(3, 7, g);

		for( core::Size ii(0); ii < 8; ++ii ) {
			TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( ii, g ), 4);
		}
	}

	void test_indole() {
		Graph g(9);
		//ring1
		boost::add_edge(0, 1, g);
		boost::add_edge(1, 2, g);
		boost::add_edge(2, 3, g);
		boost::add_edge(3, 4, g);
		boost::add_edge(4, 5, g);
		boost::add_edge(5, 0, g);
		//ring 2
		boost::add_edge(3, 6, g);
		boost::add_edge(6, 7, g);
		boost::add_edge(7, 8, g);
		boost::add_edge(8, 4, g);

		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 0, g ), 6);
		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 1, g ), 6);
		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 2, g ), 6);
		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 3, g ), 5);
		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 4, g ), 5);
		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 5, g ), 6);
		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 6, g ), 5);
		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 7, g ), 5);
		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 8, g ), 5);

		//Now test the ring limit
		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 1, g, 5), 999999);
	}
	void test_early_out() {
		//This tries to test against a possible "early out" bug. 
		//You can get n+1 rings before you get n rings, even with BFS
		Graph g(9);
		// Welcome ... windowpane!
		boost::add_edge(0, 1, g);
		boost::add_edge(0, 2, g);
		boost::add_edge(0, 3, g);
		boost::add_edge(0, 4, g);
		boost::add_edge(1, 5, g);
		boost::add_edge(2, 6, g);
		boost::add_edge(3, 7, g);
		boost::add_edge(4, 8, g);
		boost::add_edge(5, 2, g);
		boost::add_edge(6, 3, g);
		boost::add_edge(7, 4, g);
		boost::add_edge(8, 1, g);
		// now a shortcut, which should be the last bond visited (cross fingers)
		boost::add_edge(2, 3, g);
		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 0, g ), 3);
		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 2, g ), 3);
		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 3, g ), 3);
		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 6, g ), 3);
		TS_ASSERT_EQUALS(utility::graph::smallest_ring_size( 5, g ), 4);

	}
};

}  // anonymous namespace
