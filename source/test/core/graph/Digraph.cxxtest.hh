// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/graph/Digraph.cxxtest.hh
/// @brief  test suite for core::graph::Digraph.cc
/// @author Ron Jacak (ron.jacak@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/types.hh>
#include <core/graph/Digraph.hh>

// C++ headers, for debugging your tests
#include <sstream>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

using namespace core::graph;


// --------------- Test Class --------------- //

class DigraphTests : public CxxTest::TestSuite {

public:

	// Shared data elements go here.
	DigraphOP g;


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {

		g = DigraphOP( new Digraph(10) );
		g->add_edge(1,2);
		g->add_edge(2,3);
		g->add_edge(3,4);
		g->add_edge(4,5);
		g->add_edge(5,6);
		g->add_edge(6,7);
		g->add_edge(6,10);
		g->add_edge(7,2);
		g->add_edge(6,5);

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
		TS_ASSERT( g->get_edge_exists(7,2) );
		TS_ASSERT( g->get_edge_exists(6,5) );
		TS_ASSERT( !(g->get_edge_exists(2,1)) );
		TS_ASSERT( !(g->get_edge_exists(6,8)) );
		TS_ASSERT( !(g->get_edge_exists(2,7)) );

	}

	void test_drop_all_edges_for_node() {

		g->drop_all_edges_for_node(6);

		TS_ASSERT( g->get_edge_exists(1,2) );
		TS_ASSERT( g->get_edge_exists(2,3) );
		TS_ASSERT( g->get_edge_exists(3,4) );
		TS_ASSERT( g->get_edge_exists(4,5) );
		TS_ASSERT( !(g->get_edge_exists(5,6)) );
		TS_ASSERT( !(g->get_edge_exists(6,5)) );
		TS_ASSERT( !(g->get_edge_exists(6,7)) );
		TS_ASSERT( !(g->get_edge_exists(6,10)) );

	}

	void dont_test_memory_usage() {
		platform::Size const mem(g->getTotalMemoryUsage());
		//this code attempts to account for the memory variability of platform::Size across platforms
		platform::Size const num_var(74);
		platform::Size const oldsize_64(2296);
		platform::Size const oldsize_32(1296);
		platform::Size const newsize_64(oldsize_64 + num_var*(sizeof(platform::Size)-4));

		platform::Size const newsize_32(oldsize_32 + num_var*(sizeof(platform::Size)-4));

		std::cout << mem << " " << newsize_64 << " " << newsize_32 << std::endl;
		TS_ASSERT((mem == newsize_32) || (mem == newsize_64));
		//TS_ASSERT((mem == 1292) || (mem == 2288)); //32 and 64 bits, respectively //old
		//TS_ASSERT((mem == 1292) || (mem == 2584)); //32 and 64 bits, respectively //new

	}

	void test_find_edge() {

		DirectedEdge* edge_6_7 = g->find_edge(6,7);
		// How do you test that the edge you got back is correct? Non-zero?
		TS_ASSERT( edge_6_7 != 0 );
		TS_ASSERT( edge_6_7->get_tail_node_ind() == 6 );
		TS_ASSERT( edge_6_7->get_head_node_ind() == 7 );

		DirectedEdge * edge_7_6 = g->find_edge( 7, 6 );
		TS_ASSERT( ! edge_7_6 );

		g->delete_edge( edge_6_7 );
		TS_ASSERT( !(g->get_edge_exists(6,7)) );

		DirectedEdge* edge_6_8 = g->find_edge(6,8);
		TS_ASSERT( edge_6_8 == 0 );

		DirectedEdge* edge_5_6 = g->find_edge(5,6);
		TS_ASSERT( edge_5_6 != 0 );
		TS_ASSERT( edge_5_6->get_tail_node_ind() == 5 );
		TS_ASSERT( edge_5_6->get_head_node_ind() == 6 );

		DirectedEdge* edge_6_5 = g->find_edge(6,5);
		TS_ASSERT( edge_6_5 != 0 );
		TS_ASSERT( edge_6_5->get_tail_node_ind() == 6 );
		TS_ASSERT( edge_6_5->get_head_node_ind() == 5 );

	}

	void test_set_num_nodes() {

		// the set_num_nodes() call should delete the graph setup by the fixture. but test that too!
		TS_ASSERT_EQUALS(g->num_nodes(), 10);
		g->set_num_nodes(5);
		TS_ASSERT_EQUALS(g->num_nodes(), 5);

	}

	void test_operator_assignment() {

		DigraphOP g2( new Digraph() );
		*g2 = *g;

		TS_ASSERT( g2->get_edge_exists(1,2) );
		TS_ASSERT( g2->get_edge_exists(2,3) );
		TS_ASSERT( g2->get_edge_exists(3,4) );
		TS_ASSERT( g2->get_edge_exists(4,5) );
		TS_ASSERT( g2->get_edge_exists(5,6) );
		TS_ASSERT( g2->get_edge_exists(6,7) );
		TS_ASSERT( g2->get_edge_exists(6,10) );
		TS_ASSERT( g2->get_edge_exists(7,2) );
		TS_ASSERT( g2->get_edge_exists(6,5) );

		// check that changes made to g2 have no effect on g
		DirectedEdge* edge_5_6 = g2->find_edge( 5, 6 );
		g->delete_edge( edge_5_6 );
		TS_ASSERT( ! g2->get_edge_exists(5,6) );
		TS_ASSERT( g->get_edge_exists(5,6) );

		// check that changes made to g have no effect on g2
		DirectedEdge* edge_4_5 = g->find_edge( 4, 5 );
		g->delete_edge( edge_4_5 );
		TS_ASSERT( ! g->get_edge_exists(4,5) );
		TS_ASSERT( g2->get_edge_exists(4,5) );


	}

	void test_copy_constructor() {

		DigraphOP g2( new Digraph( *g ) );

		TS_ASSERT( g2->get_edge_exists(1,2) );
		TS_ASSERT( g2->get_edge_exists(2,3) );
		TS_ASSERT( g2->get_edge_exists(3,4) );
		TS_ASSERT( g2->get_edge_exists(4,5) );
		TS_ASSERT( g2->get_edge_exists(5,6) );
		TS_ASSERT( g2->get_edge_exists(6,7) );
		TS_ASSERT( g2->get_edge_exists(6,10) );
		TS_ASSERT( g2->get_edge_exists(7,2) );
		TS_ASSERT( g2->get_edge_exists(6,5) );

		// check that changes made to g2 have no effect on g
		DirectedEdge* edge_5_6 = g2->find_edge( 5, 6 );
		g2->delete_edge( edge_5_6 );
		TS_ASSERT( ! g2->get_edge_exists(5,6) );
		TS_ASSERT( g->get_edge_exists(5,6) );

		// check that changes made to g have no effect on g2
		DirectedEdge* edge_4_5 = g->find_edge( 4, 5 );
		g->delete_edge( edge_4_5 );
		TS_ASSERT( ! g->get_edge_exists(4,5) );
		TS_ASSERT( g2->get_edge_exists(4,5) );


	}


	void test_get_node() {

		DirectedNode* node6 = g->get_node(6);

		TS_ASSERT( node6->find_edge_from(5) != 0 );
		TS_ASSERT( node6->find_edge_to(7) != 0 );
		TS_ASSERT( node6->find_edge_to(10) != 0 );
	}


	void test_num_neighbors_including_self() {

		Digraph g = Digraph(1);
		TS_ASSERT_EQUALS(g.get_node(1)->num_neighbors_counting_self(), 1);
	}


	void test_serialize_graph() {
		TS_ASSERT( true );
#ifdef SERIALIZATION

		using core::Size;

		Digraph g = Digraph( 10 );
		std::set< std::pair< Size, Size > > edges;
		edges.insert( std::make_pair( 1, 2 ));
		edges.insert( std::make_pair( 2, 3 ));
		edges.insert( std::make_pair( 3, 4 ));
		edges.insert( std::make_pair( 4, 3 ));
		edges.insert( std::make_pair( 5, 6 ));
		edges.insert( std::make_pair( 6, 7 ));
		edges.insert( std::make_pair( 7, 6 ));
		edges.insert( std::make_pair( 7, 8 ));
		edges.insert( std::make_pair( 8, 9 ));
		edges.insert( std::make_pair( 9, 10 ));
		edges.insert( std::make_pair( 1, 5 ));
		edges.insert( std::make_pair( 4, 9 ));
		for ( std::set< std::pair< Size, Size > >::const_iterator iter = edges.begin(), iter_end = edges.end(); iter != iter_end; ++iter ) {
			g.add_edge( iter->first, iter->second );
		}

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( g );
		}

		std::istringstream iss( oss.str() );
		Digraph g2;
		{
			cereal::BinaryInputArchive arc( iss );
			arc( g2 );
		}

		TS_ASSERT_EQUALS( g2.num_nodes(), 10 );
		for ( Size ii = 1; ii <= 10; ++ii ) {
			for ( Size jj = ii+1; jj <= 10; ++jj ) {
				if ( std::find( edges.begin(), edges.end(), std::make_pair( ii, jj ) ) == edges.end() ) {
					TS_ASSERT( ! g2.get_edge_exists( ii, jj ) );
				} else {
					TS_ASSERT( g2.get_edge_exists( ii, jj ) );
				}
			}
		}
#endif
	}

	void test_serialize_graph_pointer() {
		TS_ASSERT( true );
#ifdef SERIALIZATION
		using core::Size;

		DigraphOP g( new  Digraph( 10 ));
		std::set< std::pair< Size, Size > > edges;
		edges.insert( std::make_pair( 1, 2 ));
		edges.insert( std::make_pair( 2, 3 ));
		edges.insert( std::make_pair( 3, 4 ));
		edges.insert( std::make_pair( 4, 5 ));
		edges.insert( std::make_pair( 5, 6 ));
		edges.insert( std::make_pair( 6, 7 ));
		edges.insert( std::make_pair( 7, 8 ));
		edges.insert( std::make_pair( 8, 9 ));
		edges.insert( std::make_pair( 9, 10 ));
		edges.insert( std::make_pair( 1, 5 ));
		edges.insert( std::make_pair( 4, 9 ));
		for ( std::set< std::pair< Size, Size > >::const_iterator iter = edges.begin(), iter_end = edges.end(); iter != iter_end; ++iter ) {
			g->add_edge( iter->first, iter->second );
		}

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( g );
		}

		std::istringstream iss( oss.str() );
		DigraphOP g2;
		{
			cereal::BinaryInputArchive arc( iss );
			arc( g2 );
		}

		TS_ASSERT_EQUALS( g2->num_nodes(), 10 );
		for ( Size ii = 1; ii <= 10; ++ii ) {
			for ( Size jj = ii+1; jj <= 10; ++jj ) {
				if ( std::find( edges.begin(), edges.end(), std::make_pair( ii, jj ) ) == edges.end() ) {
					TS_ASSERT( ! g2->get_edge_exists( ii, jj ) );
				} else {
					TS_ASSERT( g2->get_edge_exists( ii, jj ) );
				}
			}
		}
#endif
	}

	void test_topological_sort_not_a_DAG() {
		std::pair< std::list< platform::Size >, bool  > result = topological_sort( *g );
		TS_ASSERT_EQUALS( result.second, false );
	}

	void test_topological_sort_is_a_DAG() {
		using platform::Size;
		g->delete_edge( g->find_edge(7,2) );
		g->delete_edge( g->find_edge(6,5) );
		std::pair< std::list< Size >, bool  > result = topological_sort( *g );
		// 1. This is a DAG.
		TS_ASSERT_EQUALS( result.second, true );

		// 2. All nodes are listed in the topological sorting
		utility::vector1< Size > nodes_in_toposort( 10, 0 );
		for ( std::list< Size >::const_iterator iter = result.first.begin(); iter != result.first.end(); ++iter ) {
			nodes_in_toposort[ *iter ] = 1;
		}
		for ( Size ii = 1; ii <= 10; ++ii ) {
			TS_ASSERT_EQUALS( nodes_in_toposort[ ii ], 1 );
		}

		// 3. Total order is as follows, for this algorithm, though
		// there are multiple total orders that represent g's
		// partial order equally well.
		int correct[] = { 9, 8, 1, 2, 3, 4, 5, 6, 10, 7 };
		utility::vector1< Size > correct_topo_v1( 10 );
		for ( int ii = 0; ii < 10; ++ii ) { correct_topo_v1[ ii+1 ] = correct[ ii ]; }
		utility::vector1< Size > actual_topo_v1( 10 );
		std::copy( result.first.begin(), result.first.end(), actual_topo_v1.begin() );
		TS_ASSERT_EQUALS( correct_topo_v1, actual_topo_v1 );
	}

};


