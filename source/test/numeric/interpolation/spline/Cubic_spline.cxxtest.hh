// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/interpolation/spline
/// @brief  test suite for numeric::interpolation::spline::Bicubic_spline
/// @author Steven Combs (steven.combs@vanderbilt.edu)
/// This tests the functions that are in the cubic spline class except for
/// the e_periodic steps.


// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <numeric/interpolation/spline/Cubic_spline.hh>
#include <numeric/interpolation/spline/Cubic_spline.fwd.hh>
#include <numeric/MathVector_operations.hh>



// --------------- Test Class --------------- //


class Cubic_spline_tests : public CxxTest::TestSuite {



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

	void test_cubic_spline_data_access(){


		const numeric::Real values[] =
		{
				26, 3, 1, 2, 1, 3, 6, 3, 8,
				2, 7, 8, 3, 4, 2, 1, 2, 5,
				30, 0, 2, 4, 6, 3, 4, 3, 3,
				4, 11, 5, 8, 5, 2, 0, 2, 2
		};
		const numeric::MathVector<numeric::Real> input_values(36, values);

		numeric::interpolation::spline::CubicSpline naturalspline;


		naturalspline.train(numeric::interpolation::spline::e_Natural, -180, 10, input_values, std::pair<numeric::Real, numeric::Real>(10,10));


		TS_ASSERT_DELTA(26.8509, naturalspline.F(-180.30) , .001);
		TS_ASSERT_DELTA(-2.83624,  naturalspline.dF(-180.30), .001 );

		TS_ASSERT_DELTA(1.15858,  naturalspline.F(180.30), .001 );
		TS_ASSERT_DELTA(-0.0816911,  naturalspline.dF(180.30), .001 );


		TS_ASSERT_DELTA(25.1493,  naturalspline.F(-179.70), .001 );
		TS_ASSERT_DELTA(-2.83479,  naturalspline.dF(-179.70), .001 );

		TS_ASSERT_DELTA(1.2076,  naturalspline.F(179.70), .001 );
		TS_ASSERT_DELTA(-0.0816911,  naturalspline.dF(179.70), .001 );



		TS_ASSERT_EQUALS(36,  naturalspline.get_dsecox().size() );
		TS_ASSERT_EQUALS(-180,  naturalspline.get_start() );
		TS_ASSERT_EQUALS(10, naturalspline.get_delta());
		TS_ASSERT_EQUALS(36, naturalspline.get_values().size());





	}


};
