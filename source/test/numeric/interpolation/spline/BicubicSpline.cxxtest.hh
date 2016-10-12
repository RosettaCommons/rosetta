// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/interpolation/spline
/// @brief  test suite for numeric::interpolation::spline::BicubicSpline
/// @author Steven Combs (steven.combs@vanderbilt.edu)
/// This tests the functions that are in the bicubic spline class except for
/// the e_periodic steps.


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <numeric/interpolation/spline/BicubicSpline.hh>
#include <numeric/interpolation/spline/CubicSpline.hh>
#include <numeric/interpolation/spline/CubicSpline.fwd.hh>
#include <numeric/MathVector_operations.hh>


// --------------- Test Class --------------- //


class BicubicSpline_tests : public CxxTest::TestSuite {


public:
	//shared data


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp()
	{

	} //Match contents of Histogram_sample.hist


	// Shared finalization goes here.
	void tearDown() {

	}

	void test_bicubic_spline_data_access(){

		numeric::Real values[] =
		{
				26,  3,   1,  2, 1, 3, 6, 3, 8, 2, 7, 8, 3, 4, 2, 1, 2,  5,  30,  0,  2,  4, 6, 3, 4, 3, 3, 4, 11, 5, 8, 5, 2, 0, 2,  2,
				2,   0,   0,  1, 2, 0, 2, 1, 1, 3, 3, 0, 1, 0, 0, 1, 0,  0,  4,   0,  0,  0, 1, 3, 0, 5, 0, 0, 1,  1, 1, 0, 0, 3, 0,  0,
				1,   0,   0,  0, 0, 0, 1, 0, 0, 1, 2, 3, 0, 0, 0, 1, 0,  0,  76,  0,  0,  0, 0, 0, 0, 3, 2, 2, 0,  0, 0, 0, 0, 0, 0,  0,
				476, 107, 19, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 9, 3, 10, 14, 447, 14, 13, 4, 9, 5, 1, 3, 0, 0, 0,  0, 0, 0, 1, 5, 14, 104,
				1,   1,   4,  2, 1, 2, 2, 1, 1, 0, 1, 0, 2, 1, 0, 0, 0,  1,  35,  0,  0,  0, 4, 0, 0, 1, 4, 3, 0,  0, 0, 2, 2, 1, 1,  1
		};


		numeric::MathMatrix<numeric::Real> input_values(5,36, values);

		numeric::interpolation::spline::BorderFlag behavior[2] = {numeric::interpolation::spline::e_Natural, 			numeric::interpolation::spline::e_Periodic};
		const numeric::Real start[2] = {10, -180};
		const numeric::Real delta[2] = {10, 10};
		const bool lin_cont[2] ={true, true};

		const std::pair<numeric::Real, numeric::Real> first_be[2] = {std::pair<numeric::Real, numeric::Real>(10,10), std::pair<numeric::Real, numeric::Real>(10,10) };


		numeric::interpolation::spline::BicubicSpline naturalspline;

		//numeric::interpolation::spline::BicubicSpline *testspline;

		//testspline->train(behavior, start, delta, input_values, lin_cont, first_be);

		naturalspline.train(behavior, start, delta, input_values, lin_cont, first_be);


		TS_ASSERT_EQUALS(36, naturalspline.get_dsecox().get_number_cols());
		TS_ASSERT_EQUALS(36, naturalspline.get_dsecoy().get_number_cols());
		TS_ASSERT_EQUALS(36, naturalspline.get_dsecoxy().get_number_cols());


		TS_ASSERT_EQUALS(33.5, naturalspline.F(numeric::MakeVector(3.0, 0.0)));
		TS_ASSERT_EQUALS(-0.5, naturalspline.dFdx(numeric::MakeVector(3.0, 0.0)));
		TS_ASSERT_LESS_THAN(-0.391435, naturalspline.dFdy(numeric::MakeVector(3.0, 0.0)));


	}


};
