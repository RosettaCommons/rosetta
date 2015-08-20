// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/graph/DisjointSets.cxxtest.hh
/// @brief  test suite for the graph::DisjointSets class
/// @author Ron Jacak

// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/graph/Graph.hh>
#include <core/graph/DisjointSets.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

// Utility Headers

// Numeric headers

// Test headers
#include <test/core/init_util.hh>

//Auto Headers


// Auto Headers

static basic::Tracer TR("test.core.graph.disjointsets");

using namespace core;


// --------------- Test Class --------------- //

class DisjointSetsTests : public CxxTest::TestSuite {

public:

	// Shared data elements go here.
	core::graph::GraphOP g;

	// --------------- Test Fixture --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case.

	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	// All memory allocated via OPs; objects should destroy themselves so nothing else to do here.
	void tearDown() {}


public:

	// --------------- Test Cases --------------- //

	/// @details
	/// A trivial example that shows how to use various functions in the DisjointSets class
	///
	void test_disjoint_sets_example1() {
		//TR << "Running test_disjoint_sets_example1..." << std::endl;

		g = core::graph::GraphOP( new graph::Graph(5) );
		g->add_edge(1,2);
		g->add_edge(1,5);
		g->add_edge(4,5);

		// it's faster to know ahead of time how many elements there are, but not necessary
		graph::DisjointSets ds( 5 );

		// iterate over all edges
		for ( graph::Graph::EdgeListIter iter = g->edge_list_begin(); iter != g->edge_list_end(); ++iter ) {
			if ( ds.ds_find( (*iter)->get_first_node_ind() ) != ds.ds_find( (*iter)->get_second_node_ind() ) ) {
				ds.ds_union( (*iter)->get_first_node_ind(), (*iter)->get_second_node_ind() );
			}
		}

		// print out all the connected components as found by the union-find algorithm
		std::map< Size, utility::vector1< Size > >::iterator it;
		std::map< Size, utility::vector1< Size > > sets = ds.sets();

		/*TR << "sets: [ " << std::endl;
		for ( it = sets.begin() ; it != sets.end(); it++ ) {
			TR << "representative: " << (*it).first << " => nodes in connected component: [ ";
			for ( Size ii=1; ii <= (*it).second.size(); ++ii ) {
				TR << (*it).second[ ii ] << ", ";
			}
			TR << "]" << std::endl;
		}
		TR << "]" << std::endl;*/

		std::map< Size, utility::vector1< Size > > correct_sets;
		utility::vector1< Size > v;
		v.push_back( 3 );
		correct_sets[ 3 ] = v;
		v.clear();
		v.push_back( 1 ); v.push_back( 2 ); v.push_back( 4 ); v.push_back( 5 );
		correct_sets[ 2 ] = v; // apl -- r37173 changes the order in which edges are added to a graph

		TS_ASSERT_EQUALS( sets, correct_sets );

	}

};

