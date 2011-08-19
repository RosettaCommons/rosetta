// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/util.cxxtest.hh
/// @brief  test suite for numeric utility functions
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <numeric/util.hh>
#include <utility/vector1.hh>

// --------------- Test Class --------------- //

class NumericUtilTests : public CxxTest::TestSuite {

	public:

	numeric::Real delta;

	void setUp() {
		delta = 0.0001;
	}

	// Shared finalization goes here.
	void tearDown() {}


	/// @brief median test
	void test_median() {
		using utility::vector1;
		using namespace numeric;

		vector1< Real > values;
		// 1-10 in random order
		values.push_back( 1 );
		values.push_back( 3 );
		values.push_back( 6 );
		values.push_back( 4 );
		values.push_back( 5 );
		values.push_back( 2 );
		values.push_back( 10 );
		values.push_back( 8 );
		values.push_back( 7 );
		values.push_back( 9 );

		TS_ASSERT_DELTA( median( values ), 5, delta );
		values.push_back( 11 );
		TS_ASSERT_DELTA( median( values ), 5.5, delta );
	}

	void test_mean() {
		using utility::vector1;
		using namespace numeric;
		vector1< Real > values;
		// 1-10 in random order
		values.push_back( 1 );
		values.push_back( 3 );
		values.push_back( 6 );
		values.push_back( 4 );
		values.push_back( 5 );
		values.push_back( 2 );
		values.push_back( 10 );
		values.push_back( 8 );
		values.push_back( 7 );
		values.push_back( 9 );
		TS_ASSERT_DELTA( mean( values ), 5.5, delta );
		values.push_back( 11 );
		TS_ASSERT_DELTA( mean( values ), 6.0, delta );
	}

	void test_log() {
		TS_ASSERT_DELTA(numeric::log(2,2), 1, delta);
		TS_ASSERT_DELTA(numeric::log(4,2), 2, delta);
		TS_ASSERT_DELTA(numeric::log(8,2), 3, delta);

		TS_ASSERT_DELTA(numeric::log(8,8), 1, delta);
		TS_ASSERT_DELTA(numeric::log(64,8), 2, delta);
		TS_ASSERT_DELTA(numeric::log(512,8), 3, delta);
	}
}; // NumericUtilTests
