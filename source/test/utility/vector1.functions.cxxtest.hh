// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/vector1.functions.cxxtest.hh
/// @brief  test suite for the functions that live in utility/vector1.functions.hh
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>

//#include <core/init_util.hh>

// Package Headers
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>

#include <numeric/constants.hh>

#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

class Vector1FunctionsTests : public CxxTest::TestSuite
{
public:

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_nmers_of_6_choose_3() { 
		utility::vector1< int > ranges( 6 );
		ranges[ 1  ] = 1;
		ranges[ 2  ] = 20;
		ranges[ 3  ] = 23;
		ranges[ 4  ] = 54;
		ranges[ 5  ] = 55;
		ranges[ 6  ] = 70;
		auto nmers = utility::nmers_of( ranges, 3 );

		// Each one should be 3 long.
		for ( auto const & elem : nmers ) {
			TS_ASSERT_EQUALS( elem.size(), 3 );
		}
		// There should be 6 choose 3 of them, or 20
		TS_ASSERT_EQUALS( nmers.size(), 20 );
	}

	void test_binary_search_ranges_even_n_elements() {
		utility::vector1< platform::Size > ranges( 10 );
		ranges[ 1  ] = 1;
		ranges[ 2  ] = 20;
		ranges[ 3  ] = 23;
		ranges[ 4  ] = 54;
		ranges[ 5  ] = 55;
		ranges[ 6  ] = 70;
		ranges[ 7  ] = 70;
		ranges[ 8  ] = 92;
		ranges[ 9  ] = 100;
		ranges[ 10 ] = 101;

		utility::vector1< std::pair< platform::Size, platform::Size > > queries_and_answers = {
			{1, 1}, {2,1}, {12,1}, {19,1}, {20,2}, {21,2}, {22,2},
			{23,3}, {24,3}, {37,3}, {53,3}, {54,4}, {55,5}, {56,5},
			{59,5}, {68,5}, {69,5}, {70,7}, {71,7}, {91,7}, {92,8},
			{97,8}, {99,8}, {100,9}, {101,10}, {120,10}, {1150,10} };
		for ( auto q_and_a : queries_and_answers ) {
			platform::Size index = utility::binary_search_ranges( ranges, q_and_a.first );
			TS_ASSERT_EQUALS( index, q_and_a.second );
			if ( index != q_and_a.second ) {
				std::cerr << "Wrong index returned for query " << q_and_a.first << "; expected " <<
					q_and_a.second << " but got " << index << std::endl;
			}
		}

	}

	void test_binary_search_ranges_odd_n_elements() {
		utility::vector1< platform::Size > ranges( 9 );
		ranges[ 1  ] = 101;
		ranges[ 2  ] = 120;
		ranges[ 3  ] = 123;
		ranges[ 4  ] = 154;
		ranges[ 5  ] = 155;
		ranges[ 6  ] = 170;
		ranges[ 7  ] = 170;
		ranges[ 8  ] = 192;
		ranges[ 9  ] = 200;

		utility::vector1< std::pair< platform::Size, platform::Size > > queries_and_answers = {
			{101,1}, {102,1}, {112,1}, {119,1}, {120,2}, {121,2}, {122,2},
			{123,3}, {124,3}, {137,3}, {153,3}, {154,4}, {155,5}, {156,5},
			{159,5}, {168,5}, {169,5}, {170,7}, {171,7}, {191,7}, {192,8},
			{197,8}, {199,8}, {200,9}, {201,9}, {220,9}, {1150,9} };
		for ( auto q_and_a : queries_and_answers ) {
			platform::Size index = utility::binary_search_ranges( ranges, q_and_a.first );
			TS_ASSERT_EQUALS( index, q_and_a.second );
			if ( index != q_and_a.second ) {
				std::cerr << "Wrong index returned for query " << q_and_a.first << "; expected " <<
					q_and_a.second << " but got " << index << std::endl;
			}
		}

	}

};

