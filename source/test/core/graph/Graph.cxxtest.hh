// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/graph/Graph.cxxtest.hh
/// @brief  test suite for core::graph::Graph.cc
/// @author Ron Jacak (ron.jacak@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/graph/Graph.hh>

// C++ headers, for debugging your tests
// AUTO-REMOVED #include <iostream>

//Auto Headers



using namespace core::graph;


// --------------- Test Class --------------- //

class GraphTests : public CxxTest::TestSuite {

	public:

	// Shared data elements go here.
	GraphOP g;


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {

		g = GraphOP( new Graph(10) );
		g->add_edge(1,2);
		g->add_edge(2,3);
		g->add_edge(3,4);
		g->add_edge(4,5);
		g->add_edge(5,6);
		g->add_edge(6,7);
		g->add_edge(6,10);

	}

	// Shared finalization goes here.
	void tearDown() {
		// Being a smart pointer, g should be destructed and "free'd" correctly, but a g->delete_everything()
		// could be placed here, if desired.
	}


	// --------------- Test Cases --------------- //

	void test_add_edge() {

		TS_ASSERT( g->get_edge_exists(1,2) );
		TS_ASSERT( g->get_edge_exists(2,3) );
		TS_ASSERT( g->get_edge_exists(3,4) );
		TS_ASSERT( g->get_edge_exists(4,5) );
		TS_ASSERT( g->get_edge_exists(5,6) );
		TS_ASSERT( g->get_edge_exists(6,7) );
		TS_ASSERT( g->get_edge_exists(6,10) );
		TS_ASSERT( !(g->get_edge_exists(1,3)) );
		TS_ASSERT( !(g->get_edge_exists(6,8)) );

	}

	void test_drop_all_edges_for_node() {

		g->drop_all_edges_for_node(6);

		TS_ASSERT( g->get_edge_exists(1,2) );
		TS_ASSERT( g->get_edge_exists(2,3) );
		TS_ASSERT( g->get_edge_exists(3,4) );
		TS_ASSERT( g->get_edge_exists(4,5) );
		TS_ASSERT( !(g->get_edge_exists(5,6)) );
		TS_ASSERT( !(g->get_edge_exists(6,7)) );
		TS_ASSERT( !(g->get_edge_exists(6,10)) );

	}

	void test_memory_usage() {
		platform::Size const mem(g->getTotalMemoryUsage());
		//this code attempts to account for the memory variability of platform::Size across platforms
		platform::Size const num_var(74);
#ifdef PTR_MODERN
		platform::Size const oldsize_64(2296);
		platform::Size const oldsize_32(1296);
#else
		platform::Size const oldsize_64(2288);
		platform::Size const oldsize_32(1292);
#endif
		platform::Size const newsize_64(oldsize_64 + num_var*(sizeof(platform::Size)-4));

		platform::Size const newsize_32(oldsize_32 + num_var*(sizeof(platform::Size)-4));

		//std::cout << mem << " " << newsize_64 << " " << newsize_32 << std::endl;
		TS_ASSERT((mem == newsize_32) || (mem == newsize_64));
		//TS_ASSERT((mem == 1292) || (mem == 2288)); //32 and 64 bits, respectively //old
		//TS_ASSERT((mem == 1292) || (mem == 2584)); //32 and 64 bits, respectively //new

	}

	void test_find_edge() {

		Edge* edge_6_7 = g->find_edge(6,7);
		// How do you test that the edge you got back is correct? Non-zero?
		TS_ASSERT( edge_6_7 != 0 );
		g->delete_edge( edge_6_7 );
		TS_ASSERT( !(g->get_edge_exists(6,7)) );

		Edge* edge_6_8 = g->find_edge(6,8);
		TS_ASSERT( edge_6_8 == 0 );

	}

	void test_set_num_nodes() {

		// the set_num_nodes() call should delete the graph setup by the fixture. but test that too!
		TS_ASSERT_EQUALS(g->num_nodes(), 10);
		g->set_num_nodes(5);
		TS_ASSERT_EQUALS(g->num_nodes(), 5);

		platform::Size const mem(g->getTotalMemoryUsage());
		//this code attempts to account for the memory variability of platform::Size across platforms
		platform::Size const num_var(30);
#ifdef PTR_MODERN
		platform::Size const oldsize_64(656);
		platform::Size const oldsize_32(388);
#else
		platform::Size const oldsize_64(648);
		platform::Size const oldsize_32(384);
#endif
		platform::Size const newsize_64(oldsize_64 + num_var*(sizeof(platform::Size)-4));

		platform::Size const newsize_32(oldsize_32 + num_var*(sizeof(platform::Size)-4));

		//std::cout << mem << " " << newsize_64 << " " << newsize_32 << std::endl;
		TS_ASSERT((mem == newsize_32) || (mem == newsize_64));
		//TS_ASSERT((mem == 384) || (mem == 648)); //32 and 64 bits, respectively //old
		//TS_ASSERT((mem == 384) || (mem == 768)); //32 and 64 bits, respectively //new
	}

	void test_operator_assignment() {

		GraphOP g2( new Graph() );
		*g2 = *g;

		TS_ASSERT( g2->get_edge_exists(1,2) );
		TS_ASSERT( g2->get_edge_exists(2,3) );
		TS_ASSERT( g2->get_edge_exists(3,4) );
		TS_ASSERT( g2->get_edge_exists(4,5) );
		TS_ASSERT( g2->get_edge_exists(5,6) );
		TS_ASSERT( g2->get_edge_exists(6,7) );
		TS_ASSERT( g2->get_edge_exists(6,10) );

		// check that changes made to g2 have no effect on g
		Edge* edge_5_6 = g2->find_edge( 5, 6 );
		g->delete_edge( edge_5_6 );
		TS_ASSERT( ! g2->get_edge_exists(5,6) );
		TS_ASSERT( g->get_edge_exists(5,6) );

		// check that changes made to g have no effect on g2
		Edge* edge_4_5 = g->find_edge( 4, 5 );
		g->delete_edge( edge_4_5 );
		TS_ASSERT( ! g->get_edge_exists(4,5) );
		TS_ASSERT( g2->get_edge_exists(4,5) );


	}

	void test_copy_constructor() {

		GraphOP g2( new Graph( *g ) );

		TS_ASSERT( g2->get_edge_exists(1,2) );
		TS_ASSERT( g2->get_edge_exists(2,3) );
		TS_ASSERT( g2->get_edge_exists(3,4) );
		TS_ASSERT( g2->get_edge_exists(4,5) );
		TS_ASSERT( g2->get_edge_exists(5,6) );
		TS_ASSERT( g2->get_edge_exists(6,7) );
		TS_ASSERT( g2->get_edge_exists(6,10) );

		// check that changes made to g2 have no effect on g
		Edge* edge_5_6 = g2->find_edge( 5, 6 );
		g2->delete_edge( edge_5_6 );
		TS_ASSERT( ! g2->get_edge_exists(5,6) );
		TS_ASSERT( g->get_edge_exists(5,6) );

		// check that changes made to g have no effect on g2
		Edge* edge_4_5 = g->find_edge( 4, 5 );
		g->delete_edge( edge_4_5 );
		TS_ASSERT( ! g->get_edge_exists(4,5) );
		TS_ASSERT( g2->get_edge_exists(4,5) );


	}

	void test_loops() {

		// The following triggers an assertion which can't be caught and causes testing to stop.
		g->add_edge(7,7);
		TS_ASSERT( g->get_edge_exists(7,7) );
		Node* node7 = g->get_node( 7 );
		TS_ASSERT( node7->num_edges() == 3 );
		Edge* edge_7_7 = g->find_edge( 7,7 );
		g->delete_edge( edge_7_7 );

		TS_ASSERT( ! g->get_edge_exists(7,7) );
		//print_edges( g->get_node( 7 ) );
		//g->print_vertices();

	}

	void test_get_node() {

		Node* node6 = g->get_node(6);

		TS_ASSERT( node6->find_edge(5) != 0 );
		TS_ASSERT( node6->find_edge(7) != 0 );
		TS_ASSERT( node6->find_edge(10) != 0 );
	}


	void test_num_neighbors_including_self() {

		Graph g = Graph(1);
		TS_ASSERT_EQUALS(g.get_node(1)->num_neighbors_counting_self(), 1);
		TS_ASSERT_EQUALS(g.get_node(1)->num_neighbors_counting_self(), 1);
	}
};


