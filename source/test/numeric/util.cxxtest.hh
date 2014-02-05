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

class NumericUtilTests : public CxxTest::TestSuite {
 public:
	numeric::Real delta;

	void setUp() {
		delta = 0.0001;
	}

	void test_clamp() {
		TS_ASSERT_EQUALS(5, numeric::clamp<int>(10, 1, 5));  // clamp 10 to interval [1, 5]
		TS_ASSERT_EQUALS(1, numeric::clamp<int>(-1, 1, 5));  // clamp -1 to interval [1, 5]

    TS_ASSERT_EQUALS(5.0, numeric::clamp<double>(10.0, 1.0, 5.0));  // clamp 10.0 to interval [1.0, 5.0]
		TS_ASSERT_EQUALS(1.0, numeric::clamp<double>(-1.0, 1.0, 5.0));  // clamp -1.0 to interval [1.0, 5.0]
	}

	void test_median() {
		using utility::vector1;
		using namespace numeric;

		vector1< Real > values;
		// 1-10 in random order
		values.push_back( 1 );
		TS_ASSERT_DELTA( median( values ), 1, delta );
		values.push_back( 3 );
		TS_ASSERT_DELTA( median( values ), 2, delta );
		values.push_back( 6 );
		TS_ASSERT_DELTA( median( values ), 3, delta );
		values.push_back( 4 );
		values.push_back( 5 );
		values.push_back( 2 );
		values.push_back( 10 );
		values.push_back( 8 );
		values.push_back( 7 );
		values.push_back( 9 );

		TS_ASSERT_DELTA( median( values ), 5.5, delta );
		values.push_back( 11 );
		TS_ASSERT_DELTA( median( values ), 6, delta );
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

  void test_min() {
    using utility::vector1;
    using namespace numeric;
    vector1< Real > values;
    // 1-10 in random order
    values.push_back( 7 );
    TS_ASSERT_DELTA( min( values ), 7, delta );
    values.push_back( 3 );
    TS_ASSERT_DELTA( min( values ), 3, delta );
    values.push_back( 6 );
    TS_ASSERT_DELTA( min( values ), 3, delta );
    values.push_back( 4 );
    values.push_back( 5 );
    values.push_back( 2 );
    values.push_back( 10 );
    values.push_back( 8 );
    values.push_back( 1 );
    values.push_back( 9 );
    TS_ASSERT_DELTA( min( values ), 1, delta );
    values.push_back( 11 );
    TS_ASSERT_DELTA( min( values ), 1, delta );
  }

  void test_max() {
    using utility::vector1;
    using namespace numeric;
    vector1< Real > values;
    // 1-10 in random order
    values.push_back( 7 );
    TS_ASSERT_DELTA( max( values ), 7, delta );
    values.push_back( 3 );
    TS_ASSERT_DELTA( max( values ), 7, delta );
    values.push_back( 6 );
    TS_ASSERT_DELTA( max( values ), 7, delta );
    values.push_back( 4 );
    values.push_back( 5 );
    values.push_back( 2 );
    values.push_back( 10 );
    values.push_back( 8 );
    values.push_back( 1 );
    values.push_back( 9 );
    TS_ASSERT_DELTA( max( values ), 10, delta );
    values.push_back( 11 );
    TS_ASSERT_DELTA( max( values ), 11, delta );
  }

	void test_log() {
		TS_ASSERT_DELTA(numeric::log(2,2), 1, delta);
		TS_ASSERT_DELTA(numeric::log(4,2), 2, delta);
		TS_ASSERT_DELTA(numeric::log(8,2), 3, delta);

		TS_ASSERT_DELTA(numeric::log(8,8), 1, delta);
		TS_ASSERT_DELTA(numeric::log(64,8), 2, delta);
		TS_ASSERT_DELTA(numeric::log(512,8), 3, delta);
	}

	void test_find_nearest() {
		using utility::vector1;
		using namespace numeric;
		vector1< Real > values;
		values.push_back(1);
		values.push_back(2);
		values.push_back(5);
		values.push_back(3);

		Real nearest = find_nearest_value<Real>(values,2.0);
		TS_ASSERT(nearest == 2.0);

		nearest = find_nearest_value(values,5.2);
		TS_ASSERT(nearest == 5.0);


	}

	void test_undefined() {
		TS_ASSERT( numeric::is_undefined( numeric::get_undefined_size() ) );
		TS_ASSERT( numeric::is_undefined( numeric::get_undefined_real() ) );
		TS_ASSERT( ! numeric::is_undefined( numeric::Size(0) ) );
		TS_ASSERT( ! numeric::is_undefined( 0.0 ) );
	}

}; // NumericUtilTests
