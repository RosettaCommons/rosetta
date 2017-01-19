// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/DiscreteIntervalEncodingTree.cxxtest.hh
/// @brief  test suite for numeric::DiscreteIntervalEncodingTree.cxxtest.hh
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <numeric/DiscreteIntervalEncodingTree.hh>

#include <utility/vector1.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("numeric.DiscreteIntervalEncodingTree_cxxtest");

typedef numeric::DiscreteIntervalEncodingTree< platform::Size > Diet;

using namespace numeric;

// --------------- Test Class --------------- //

class DiscreteIntervalEncodingTreeTests : public CxxTest::TestSuite {
public:
	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
	/// @brief Test identity transform constructor
	void test_DiscreteIntervalEncodingTree_default_constructor() {
		Diet diet;
		TS_ASSERT_EQUALS( diet.size(), 0 );
		TS_ASSERT( diet.correct() );
	}

	void test_DiscreteIntervalEncodingTree_insert_one_element() {
		Diet diet;
		diet.insert( 5 );
		TS_ASSERT_EQUALS( diet.size(), 1 );
		TS_ASSERT( diet.correct() );
	}

	void test_DiscreteIntervalEncodingTree_insert_two_elements_to_two_nodes_left() {
		Diet diet;
		diet.insert( 5 );
		diet.insert( 3 );
		TS_ASSERT_EQUALS( diet.size(), 2 );
		TS_ASSERT( diet.correct() );
	}

	void test_DiscreteIntervalEncodingTree_insert_two_elements_to_two_nodes_right() {
		Diet diet;
		diet.insert( 5 );
		diet.insert( 7 );
		TS_ASSERT_EQUALS( diet.size(), 2 );
		TS_ASSERT( diet.correct() );
	}

	void test_DiscreteIntervalEncodingTree_insert_two_elements_to_one_node_left() {
		Diet diet;
		diet.insert( 5 );
		diet.insert( 4 );
		TS_ASSERT_EQUALS( diet.size(), 1 );
		TS_ASSERT( diet.correct() );
	}

	void test_DiscreteIntervalEncodingTree_insert_two_elements_to_one_node_right() {
		Diet diet;
		diet.insert( 5 );
		diet.insert( 6 );
		TS_ASSERT_EQUALS( diet.size(), 1 );
		TS_ASSERT( diet.correct() );
	}

	void test_DiscreteIntervalEncodingTree_insert_three_elements_to_one_node_left() {
		Diet diet;
		diet.insert( 5 );
		diet.insert( 4 );
		diet.insert( 3 );
		TS_ASSERT_EQUALS( diet.size(), 1 );
		TS_ASSERT( diet.correct() );
	}

	void test_DiscreteIntervalEncodingTree_insert_three_elements_to_one_node_right() {
		Diet diet;
		diet.insert( 5 );
		diet.insert( 6 );
		diet.insert( 7 );
		TS_ASSERT_EQUALS( diet.size(), 1 );
		TS_ASSERT( diet.correct() );
	}

	void test_DiscreteIntervalEncodingTree_find_elements() {
		Diet diet;
		diet.insert( 5 );
		diet.insert( 6 );
		diet.insert( 7 );
		diet.insert( 11 );
		diet.insert( 9 );
		diet.insert( 2 );
		diet.insert( 3 );

		TS_ASSERT( diet.correct() );

		TS_ASSERT( ! diet.member( 1 )  );
		TS_ASSERT(   diet.member( 2 )  );
		TS_ASSERT(   diet.member( 3 )  );
		TS_ASSERT( ! diet.member( 4 )  );
		TS_ASSERT(   diet.member( 5 )  );
		TS_ASSERT(   diet.member( 6 )  );
		TS_ASSERT(   diet.member( 7 )  );
		TS_ASSERT( ! diet.member( 8 )  );
		TS_ASSERT(   diet.member( 9 )  );
		TS_ASSERT( ! diet.member( 10 ) );
		TS_ASSERT(   diet.member( 11 ) );
		TS_ASSERT( ! diet.member( 12 ) );
	}

	void test_DiscreteIntervalEncodingTree_inorder_ranges() {
		Diet diet;

		diet.insert( 5 );
		diet.insert( 6 );
		diet.insert( 7 );
		diet.insert( 12 );
		diet.insert( 11 );
		diet.insert( 9 );
		diet.insert( 2 );
		diet.insert( 3 );

		TS_ASSERT( diet.correct() );

		Diet::RangeList range_list = diet.ranges();
		typedef utility::vector1< std::pair< platform::Size, platform::Size > > RangeVector;
		RangeVector range_vector( range_list.size() );
		std::copy( range_list.begin(), range_list.end(), range_vector.begin() );

		TS_ASSERT_EQUALS( range_vector.size(), 4 );

		TS_ASSERT_EQUALS( range_vector[ 1 ].first,  2 );
		TS_ASSERT_EQUALS( range_vector[ 1 ].second, 3 );

		TS_ASSERT_EQUALS( range_vector[ 2 ].first,  5 );
		TS_ASSERT_EQUALS( range_vector[ 2 ].second, 7 );

		TS_ASSERT_EQUALS( range_vector[ 3 ].first,  9 );
		TS_ASSERT_EQUALS( range_vector[ 3 ].second, 9 );

		TS_ASSERT_EQUALS( range_vector[ 4 ].first,  11 );
		TS_ASSERT_EQUALS( range_vector[ 4 ].second, 12 );
	}

	void test_diet_with_insertion_order_from_jd3() {
		Diet diet;
		std::list< platform::Size > insert_order = {
			2, 3, 4, 6, 7, 8, 9,
			11, 13, 14, 15, 16, 17, 18, 19, 20,	21, 22, 24, 25, 26,
			27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
			10, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
			51, 52, 53, 54, 56, 57, 58, 59, 60,
			12 };
		for ( auto ii : insert_order ) {
			diet.insert( ii );
		}

		Diet::RangeList range_list = diet.ranges();
		for ( auto range : range_list ) {
			TR << "range: " << range.first << ", " << range.second << std::endl;
		}

		TS_ASSERT( diet.member( 12 ) );
	}

};


