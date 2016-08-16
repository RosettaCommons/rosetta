// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/sphericalVector.cxxtest.hh
/// @brief  test suite for numeric::sphericalVector
/// @author Sam DeLuca


// Test headers
#include <cxxtest/TestSuite.h>

// Package Headers
#include <numeric/sphericalVector.hh>
//#include <numeric/xyzVector.io.hh>


// --------------- Test Class --------------- //

class SphericalVectorTests : public CxxTest::TestSuite {

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
	/// @brief test construction from 1D array
	void test_sphericalVector_Construct() {

		numeric::sphericalVector<float> v( 15.0,16.0,17.0 );
		TS_ASSERT_EQUALS( v.phi(), 15.0f );
		TS_ASSERT_EQUALS( v.theta(), 16.0f );
		TS_ASSERT_EQUALS( v.radius(), 17.0f );
	}

	void test_sphericalVector_CopyConstruct() {
		numeric::sphericalVector<float> v(15.0,16.0,17.0);
		numeric::sphericalVector<float> copy(v);
		TS_ASSERT_EQUALS(v.phi(),copy.phi());
		TS_ASSERT_EQUALS(v.theta(),copy.theta());
		TS_ASSERT_EQUALS(v.radius(),copy.radius());
	}

	void test_sphericalVector_ArithmaticOperators(){
		numeric::sphericalVector<float> v( 1.0,2.0,3.0 );
		numeric::sphericalVector<float> x( 4.0,5.0,6.0 );

		numeric::sphericalVector<float> y(v+x);
		TS_ASSERT_EQUALS(y.phi(),5.0);
		TS_ASSERT_EQUALS(y.theta(),7.0);
		TS_ASSERT_EQUALS(y.radius(),9.0);

		y = v-x;
		TS_ASSERT_EQUALS(y.phi(),-3.0);
		TS_ASSERT_EQUALS(y.theta(),-3.0);
		TS_ASSERT_EQUALS(y.radius(),-3.0);

		y = v*0.5;
		TS_ASSERT_EQUALS(y.phi(),0.5);
		TS_ASSERT_EQUALS(y.theta(),1.0);
		TS_ASSERT_EQUALS(y.radius(),1.5);

		y= v/2;
		TS_ASSERT_EQUALS(y.phi(),0.5);
		TS_ASSERT_EQUALS(y.theta(),1.0);
		TS_ASSERT_EQUALS(y.radius(),1.5);
	}


};

