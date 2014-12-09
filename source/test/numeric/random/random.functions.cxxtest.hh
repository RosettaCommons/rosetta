// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/random/random.functions.cxxtest.hh
/// @brief  test suite for the functions that live in numeric/random/random.functions.hh
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>

//#include <core/init_util.hh>

// Package Headers
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>

#include <numeric/constants.hh>

#include <utility/vector1.hh>

class random_functionsTests : public CxxTest::TestSuite
{
public:

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_binary_search_cdf() {
		utility::vector1< float > cdf( 10 );
		cdf[ 1  ] = 0.000f;
		cdf[ 2  ] = 0.125f;
		cdf[ 3  ] = 0.300f;
		cdf[ 4  ] = 0.328f;
		cdf[ 5  ] = 0.331f;
		cdf[ 6  ] = 0.390f;
		cdf[ 7  ] = 0.410f;
		cdf[ 8  ] = 0.651f;
		cdf[ 9  ] = 0.799f;
		cdf[ 10 ] = 0.881f;

		platform::Size test_find_0p450 = numeric::random::binary_search_cdf( cdf, 0.450f );
		TS_ASSERT_EQUALS( test_find_0p450, 7 );

		platform::Size test_find_0p070 = numeric::random::binary_search_cdf( cdf, 0.070f );
		TS_ASSERT_EQUALS( test_find_0p070, 1 );

		platform::Size test_find_0p970 = numeric::random::binary_search_cdf( cdf, 0.970f );
		TS_ASSERT_EQUALS( test_find_0p970, 10 );

		platform::Size test_find_0p125 = numeric::random::binary_search_cdf( cdf, 0.125f );
		TS_ASSERT_EQUALS( test_find_0p125, 2 );

	}

	void test_binary_search_cdf_nentires_odd() {
		utility::vector1< float > cdf( 11 );
		cdf[ 1  ] = 0.000f;
		cdf[ 2  ] = 0.125f;
		cdf[ 3  ] = 0.300f;
		cdf[ 4  ] = 0.328f;
		cdf[ 5  ] = 0.331f;
		cdf[ 6  ] = 0.390f;
		cdf[ 7  ] = 0.410f;
		cdf[ 8  ] = 0.651f;
		cdf[ 9  ] = 0.799f;
		cdf[ 10 ] = 0.881f;
		cdf[ 11 ] = 0.952f;

		platform::Size test_find_0p450 = numeric::random::binary_search_cdf( cdf, 0.450f );
		TS_ASSERT_EQUALS( test_find_0p450, 7 );

		platform::Size test_find_0p070 = numeric::random::binary_search_cdf( cdf, 0.070f );
		TS_ASSERT_EQUALS( test_find_0p070, 1 );

		platform::Size test_find_0p90 = numeric::random::binary_search_cdf( cdf, 0.90f );
		TS_ASSERT_EQUALS( test_find_0p90, 10 );

		platform::Size test_find_0p970 = numeric::random::binary_search_cdf( cdf, 0.970f );
		TS_ASSERT_EQUALS( test_find_0p970, 11 );

		platform::Size test_find_0p125 = numeric::random::binary_search_cdf( cdf, 0.125f );
		TS_ASSERT_EQUALS( test_find_0p125, 2 );

	}

	void test_binary_search_cdf_w_zero_prob_values() {
		utility::vector1< float > cdf( 10 );
		cdf[ 1  ] = 0.000f;
		cdf[ 2  ] = 0.000f;
		cdf[ 3  ] = 0.325f;
		cdf[ 4  ] = 0.325f;
		cdf[ 5  ] = 0.450f;
		cdf[ 6  ] = 0.550f;
		cdf[ 7  ] = 0.600f;
		cdf[ 8  ] = 0.651f;
		cdf[ 9  ] = 0.651f;
		cdf[ 10 ] = 0.881f;

		platform::Size test_find_0p01 = numeric::random::binary_search_cdf( cdf, 0.01f );
		TS_ASSERT_EQUALS( test_find_0p01, 2 );

		platform::Size test_find_0p34 = numeric::random::binary_search_cdf( cdf, 0.34f );
		TS_ASSERT_EQUALS( test_find_0p34, 4 );

		platform::Size test_find_0p5 = numeric::random::binary_search_cdf( cdf, 0.5f );
		TS_ASSERT_EQUALS( test_find_0p5, 5 );

		platform::Size test_find_0p78 = numeric::random::binary_search_cdf( cdf, 0.78f );
		TS_ASSERT_EQUALS( test_find_0p78, 9 );

	}

};

