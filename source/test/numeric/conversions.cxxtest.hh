// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/conversions.cxxtest.hh
/// @brief  numeric::conversions test suite
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <numeric/conversions.hh>


// --------------- Test Class --------------- //

class ConversionsTests : public CxxTest::TestSuite {

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


	/// @brief general conversion tests
	void test_conversions_general() {

		using numeric::conversions::degrees;
		using numeric::conversions::radians;
		using numeric::conversions::to_degrees;
		using numeric::conversions::to_radians;

		float rad_f = 1.0f;
		float deg_f = 1.0f;

		TS_ASSERT_DELTA( degrees( rad_f ), 57.2957795131f, delta_percent );
		TS_ASSERT_DELTA( radians( deg_f ),  0.0174532925f, delta_percent );

		to_degrees( rad_f );
		TS_ASSERT_DELTA( rad_f, 57.2957795131f, delta_percent );
		to_radians( deg_f );
		TS_ASSERT_DELTA( deg_f,  0.0174532925f, delta_percent );

		double rad_d = 1.0;
		double deg_d = 1.0;

		to_degrees( rad_d );
		TS_ASSERT_DELTA( rad_d, 57.29577951308232087680, delta_percent );
		to_radians( deg_d );
		TS_ASSERT_DELTA( deg_d,  0.01745329251994329577, delta_percent );
	}

};


