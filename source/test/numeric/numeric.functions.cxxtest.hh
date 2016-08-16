// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/numeric.functions.cxxtest.hh
/// @brief  test suite for numericnumeric_functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <numeric/numeric.functions.hh>


// --------------- Test Class --------------- //

class NumericFunctionsTests : public CxxTest::TestSuite {

	public:

	// Shared data elements go here.
	double delta_percent;         // percentage difference for floating-point comparisons in TS_ASSERT_DELTA

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		delta_percent = 0.0001;
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
	/// @brief nearest tests
	void test_numeric_functions_nearest() {

		using numeric::nearest;
		TS_ASSERT_EQUALS( nearest< int >( 53.3 ), 53 );
		TS_ASSERT_EQUALS( nearest< int >( 53.6 ), 54 );
		TS_ASSERT_EQUALS( nearest< long >( 53.3 ), 53l );
		TS_ASSERT_EQUALS( nearest< long >( 53.6 ), 54l );
		TS_ASSERT_DELTA( nearest< float >( 53.3 ), 53.3f, delta_percent );
		TS_ASSERT_DELTA( nearest< float >( 53.6 ), 53.6f, delta_percent );
	}

	/// @brief nint tests
	void test_numeric_functions_nint() {

		using numeric::nint;
		TS_ASSERT_EQUALS( nint( 53.3 ), 53 );
		TS_ASSERT_EQUALS( nint( 53.6 ), 54 );
	}

	/// @brief mod tests
	void test_numeric_functions_mod() {

		using numeric::mod;
		TS_ASSERT_EQUALS( mod( 7, 2 ), 1 );
		TS_ASSERT_EQUALS( mod( 29, 6 ), 5 );
		TS_ASSERT_DELTA( mod( 29.0, 6.0 ), 5.0, delta_percent );
		TS_ASSERT_DELTA( mod( 29.0, 6.5 ), 3.0, delta_percent );
	}

	/// @brief remainder tests
	void test_numeric_functions_remainder() {

		// On GCC and other compilers providing C99 math functions these will call math.h versions
		// (Why check these?  Just want to check our own versions.  Remove for now.  kph, 2006-05-24)
		//	TS_ASSERT_EQUALS( remainder( 7, 2 ), -1 );
		//	TS_ASSERT_EQUALS( remainder( 29, 6 ), -1 );
		//	TS_ASSERT_EQUALS( remainder( 27, 6 ), 3 ); // Special case: n - (x/y) == .5
		//	TS_ASSERT_EQUALS( remainder( 33, 6 ), -3 ); // Special case: n - (x/y) == -.5
		//	TS_ASSERT_DELTA( remainder( 29.0, 6.0 ), -1.0, delta_percent );
		//	TS_ASSERT_DELTA( remainder( 29.0, 6.5 ), 3.0, delta_percent );
		//	TS_ASSERT_DELTA( remainder( 27.0, 6.0 ), 3.0, delta_percent ); // Special case: n - (x/y) == .5
		//	TS_ASSERT_DELTA( remainder( 33.0, 6.0 ), -3.0, delta_percent ); // Special case: n - (x/y) == -.5

		// Make sure we call the numeric_functions version
		TS_ASSERT_EQUALS( numeric::remainder( 7, 2 ), -1 );
		TS_ASSERT_EQUALS( numeric::remainder( 29, 6 ), -1 );
		TS_ASSERT_EQUALS( numeric::remainder( 27, 6 ), 3 ); // Special case: n - (x/y) == .5
		TS_ASSERT_EQUALS( numeric::remainder( 33, 6 ), -3 ); // Special case: n - (x/y) == -.5
		TS_ASSERT_DELTA( numeric::remainder( 29.0, 6.0 ), -1.0, delta_percent );
		TS_ASSERT_DELTA( numeric::remainder( 29.0, 6.5 ), 3.0, delta_percent );
		TS_ASSERT_DELTA( numeric::remainder( 27.0, 6.0 ), 3.0, delta_percent ); // Special case: n - (x/y) == .5
		TS_ASSERT_DELTA( numeric::remainder( 33.0, 6.0 ), -3.0, delta_percent ); // Special case: n - (x/y) == -.5
	}

	/// @brief gcd tests
	void test_numeric_functions_gcd() {

		using numeric::gcd;
		TS_ASSERT_EQUALS( gcd( 7, 2 ), 1 );
		TS_ASSERT_EQUALS( gcd( 28, 6 ), 2 );
		TS_ASSERT_EQUALS( gcd( 16, 12 ), 4 );
		TS_ASSERT_DELTA( gcd( 7.0, 2.0 ), 1.0, delta_percent );
		TS_ASSERT_DELTA( gcd( 28.0, 6.0 ), 2.0, delta_percent );
		TS_ASSERT_DELTA( gcd( 16.0, 12.0 ), 4.0, delta_percent );
	}

};


