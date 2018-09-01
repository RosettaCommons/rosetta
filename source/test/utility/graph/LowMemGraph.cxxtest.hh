// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/graph/Graph.cxxtest.hh
/// @brief  test suite for utility::graph::Graph.cc
/// @author Brian Coventry (ron.jacak@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/types.hh>
#include <utility/graph/LowMemGraph.hh>

// C++ headers, for debugging your tests
#include <sstream>

using namespace core;
using namespace utility::graph;


// --------------- Test Class --------------- //

class LowMemGraphTests : public CxxTest::TestSuite {

public:

	// Shared data elements go here.
	DefaultLowMemGraphOP g;


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {

		g = DefaultLowMemGraphOP( new DefaultLowMemGraph(10) );
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
		g = 0;
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

	template< class Iter >
	void
	test_iter( Iter const & begin, Iter const & end, utility::vector1<std::pair<Size,Size>> edges ) {

		// Forwards
		{
			utility::vector1<bool> found_edge( edges.size(), false );

			for ( Iter iter = begin; iter != end; ++iter ) {
				LowMemEdge const * edge = *iter;
				std::pair<Size,Size> search_pair ( edge->get_first_node_ind(), edge->get_second_node_ind() );

				Size found_index = 0;
				for ( Size ii = 1; ii <= edges.size(); ii++ ) {
					if ( edges[ii] == search_pair ) {
						TS_ASSERT( found_index == 0 );
						found_index = ii;
					}
				}

				TS_ASSERT( found_index != 0 );
				TS_ASSERT( !found_edge[ found_index ] );
				found_edge[ found_index ] = true;

			}

			for ( Size ii = 1; ii <= edges.size(); ii++ ) {
				TS_ASSERT( found_edge[ ii ] );
			}
		}
		// Backwards
		{
			utility::vector1<bool> found_edge( edges.size(), false );

			Iter iter = end;
			while ( iter != begin ) {
				--iter;

				LowMemEdge const * edge = *iter;
				std::pair<Size,Size> search_pair ( edge->get_first_node_ind(), edge->get_second_node_ind() );

				Size found_index = 0;
				for ( Size ii = 1; ii <= edges.size(); ii++ ) {
					if ( edges[ii] == search_pair ) {
						TS_ASSERT( found_index == 0 );
						found_index = ii;
					}
				}

				TS_ASSERT( found_index != 0 );
				TS_ASSERT( !found_edge[ found_index ] );
				found_edge[ found_index ] = true;

			}

			for ( Size ii = 1; ii <= edges.size(); ii++ ) {
				TS_ASSERT( found_edge[ ii ] );
			}
		}


	}

	void test_iterators() {
		utility::vector1<std::pair<Size,Size>> all_edges {
			{1, 2},
			{2, 3},
			{3, 4},
			{4, 5},
			{5, 6},
			{6, 7},
			{6, 10}
			};
		utility::vector1<std::pair<Size,Size>> five_edges {
			{4, 5},
			{5, 6}
			};
		utility::vector1<std::pair<Size,Size>> no_edges {
			};

		test_iter( g->edge_list_begin(), g->edge_list_end(), all_edges );
		test_iter( g->const_edge_list_begin(), g->const_edge_list_end(), all_edges );

		LowMemNode * node5 = g->get_node( 5 );

		test_iter( node5->edge_list_begin(*g), node5->edge_list_end(*g), five_edges );
		test_iter( node5->const_edge_list_begin(*g), node5->const_edge_list_end(*g), five_edges );

		LowMemNode * node8 = g->get_node( 8 );

		test_iter( node8->edge_list_begin(*g), node8->edge_list_end(*g), no_edges );
		test_iter( node8->const_edge_list_begin(*g), node8->const_edge_list_end(*g), no_edges );

		g->drop_all_edges();
		test_iter( g->edge_list_begin(), g->edge_list_end(), no_edges );
		test_iter( g->const_edge_list_begin(), g->const_edge_list_end(), no_edges );
	}

	void test_iterators_w_deleted_data() {
		g->drop_all_edges_for_node(1);
		g->drop_all_edges_for_node(6);


		utility::vector1<std::pair<Size,Size>> all_edges {
			{2, 3},
			{3, 4},
			{4, 5},
			};
		utility::vector1<std::pair<Size,Size>> five_edges {
			{4, 5}
			};

		test_iter( g->edge_list_begin(), g->edge_list_end(), all_edges );
		test_iter( g->const_edge_list_begin(), g->const_edge_list_end(), all_edges );

		LowMemNode * node5 = g->get_node( 5 );

		test_iter( node5->edge_list_begin(*g), node5->edge_list_end(*g), five_edges );
		test_iter( node5->const_edge_list_begin(*g), node5->const_edge_list_end(*g), five_edges );
	}

	// void test_memory_usage() {
	//  platform::Size const mem(g->getTotalMemoryUsage());
	//  //this code attempts to account for the memory variability of platform::Size across platforms
	//  platform::Size const num_var(74);
	//  platform::Size const oldsize_64(2296);
	//  platform::Size const oldsize_32(1296);
	//  platform::Size const newsize_64(oldsize_64 + num_var*(sizeof(platform::Size)-4));

	//  platform::Size const newsize_32(oldsize_32 + num_var*(sizeof(platform::Size)-4));

	//  //std::cout << mem << " " << newsize_64 << " " << newsize_32 << std::endl;
	//  TS_ASSERT((mem == newsize_32) || (mem == newsize_64));
	//  //TS_ASSERT((mem == 1292) || (mem == 2288)); //32 and 64 bits, respectively //old
	//  //TS_ASSERT((mem == 1292) || (mem == 2584)); //32 and 64 bits, respectively //new

	// }

	void test_find_edge() {

		LowMemEdge* edge_6_7 = g->find_edge(6,7);
		// How do you test that the edge you got back is correct? Non-zero?
		TS_ASSERT( edge_6_7 != 0 );
		g->delete_edge( edge_6_7 );
		TS_ASSERT( !(g->get_edge_exists(6,7)) );

		LowMemEdge* edge_6_8 = g->find_edge(6,8);
		TS_ASSERT( edge_6_8 == 0 );

	}

	// void test_set_num_nodes() {

	//  // the set_num_nodes() call should delete the graph setup by the fixture. but test that too!
	//  TS_ASSERT_EQUALS(g->num_nodes(), 10);
	//  g->set_num_nodes(5);
	//  TS_ASSERT_EQUALS(g->num_nodes(), 5);

	//  platform::Size const mem(g->getTotalMemoryUsage());
	//  //this code attempts to account for the memory variability of platform::Size across platforms
	//  platform::Size const num_var(30);
	//  platform::Size const oldsize_64(656);
	//  platform::Size const oldsize_32(388);
	//  platform::Size const newsize_64(oldsize_64 + num_var*(sizeof(platform::Size)-4));

	//  platform::Size const newsize_32(oldsize_32 + num_var*(sizeof(platform::Size)-4));

	//  //std::cout << mem << " " << newsize_64 << " " << newsize_32 << std::endl;
	//  TS_ASSERT((mem == newsize_32) || (mem == newsize_64));
	//  //TS_ASSERT((mem == 384) || (mem == 648)); //32 and 64 bits, respectively //old
	//  //TS_ASSERT((mem == 384) || (mem == 768)); //32 and 64 bits, respectively //new
	// }

	void test_operator_assignment() {

		DefaultLowMemGraphOP g2( new DefaultLowMemGraph() );
		*g2 = *g;

		TS_ASSERT( g2->get_edge_exists(1,2) );
		TS_ASSERT( g2->get_edge_exists(2,3) );
		TS_ASSERT( g2->get_edge_exists(3,4) );
		TS_ASSERT( g2->get_edge_exists(4,5) );
		TS_ASSERT( g2->get_edge_exists(5,6) );
		TS_ASSERT( g2->get_edge_exists(6,7) );
		TS_ASSERT( g2->get_edge_exists(6,10) );

		// check that changes made to g2 have no effect on g
		LowMemEdge* edge_5_6 = g2->find_edge( 5, 6 );
		g2->delete_edge( edge_5_6 );
		TS_ASSERT( ! g2->get_edge_exists(5,6) );
		TS_ASSERT( g->get_edge_exists(5,6) );

		// check that changes made to g have no effect on g2
		LowMemEdge* edge_4_5 = g->find_edge( 4, 5 );
		g->delete_edge( edge_4_5 );
		TS_ASSERT( ! g->get_edge_exists(4,5) );
		TS_ASSERT( g2->get_edge_exists(4,5) );



	}

	void test_copy_constructor() {

		DefaultLowMemGraphOP g2( new DefaultLowMemGraph( *g ) );

		TS_ASSERT( g2->get_edge_exists(1,2) );
		TS_ASSERT( g2->get_edge_exists(2,3) );
		TS_ASSERT( g2->get_edge_exists(3,4) );
		TS_ASSERT( g2->get_edge_exists(4,5) );
		TS_ASSERT( g2->get_edge_exists(5,6) );
		TS_ASSERT( g2->get_edge_exists(6,7) );
		TS_ASSERT( g2->get_edge_exists(6,10) );

		// check that changes made to g2 have no effect on g
		LowMemEdge* edge_5_6 = g2->find_edge( 5, 6 );
		g2->delete_edge( edge_5_6 );
		TS_ASSERT( ! g2->get_edge_exists(5,6) );
		TS_ASSERT( g->get_edge_exists(5,6) );

		// check that changes made to g have no effect on g2
		LowMemEdge* edge_4_5 = g->find_edge( 4, 5 );
		g->delete_edge( edge_4_5 );
		TS_ASSERT( ! g->get_edge_exists(4,5) );
		TS_ASSERT( g2->get_edge_exists(4,5) );


	}

	// void test_loops() {

	//  // The following triggers an assertion which can't be caught and causes testing to stop.
	//  g->add_edge(7,7);
	//  TS_ASSERT( g->get_edge_exists(7,7) );
	//  Node* node7 = g->get_node( 7 );
	//  TS_ASSERT( node7->num_edges() == 3 );
	//  Edge* edge_7_7 = g->find_edge( 7,7 );
	//  g->delete_edge( edge_7_7 );

	//  TS_ASSERT( ! g->get_edge_exists(7,7) );
	//  //print_edges( g->get_node( 7 ) );
	//  //g->print_vertices();

	// }

	void test_get_node() {

		LowMemNode* node6 = g->get_node(6);

		TS_ASSERT( node6->find_edge(5, *g) != 0 );
		TS_ASSERT( node6->find_edge(7, *g) != 0 );
		TS_ASSERT( node6->find_edge(10, *g) != 0 );
	}

	void test_deletion_space_savings() {

		TS_ASSERT( g->internal_edge_list_size() == 7 );

		g->add_edge( 9, 10 );

		TS_ASSERT( g->internal_edge_list_size() == 8 );

		LowMemEdge* edge_5_6 = g->find_edge( 5, 6 );
		g->delete_edge( edge_5_6 );

		TS_ASSERT( g->internal_edge_list_size() == 8 );

		g->add_edge( 8, 10 );

		TS_ASSERT( g->internal_edge_list_size() == 8 );

		utility::vector1<std::pair<Size,Size>> these_edges {
			{1, 2},
			{2, 3},
			{3, 4},
			{4, 5},
			{6, 7},
			{6, 10},
			{8, 10},
			{9, 10}
			};

		// Make sure the edge actually works
		test_iter( g->edge_list_begin(), g->edge_list_end(), these_edges );
	}


	// void test_num_neighbors_including_self() {

	//  Graph g = Graph(1);
	//  TS_ASSERT_EQUALS(g.get_node(1)->num_neighbors_counting_self(), 1);
	//  TS_ASSERT_EQUALS(g.get_node(1)->num_neighbors_counting_self(), 1);
	// }



};


