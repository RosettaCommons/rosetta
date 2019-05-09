// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  numeric/linear_algebra/minimum_bounding_ellipse.cxxtest.hh
/// @brief  Given a set of points, calculate the minimum bounding ellipse
/// @author Rebecca Alford (rfalford12@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <numeric/MathMatrix.hh>
#include <numeric/MathMatrix_operations.hh>
#include <numeric/xyzVector.hh>
#include <numeric/linear_algebra/minimum_bounding_ellipse.hh>
#include <numeric/linear_algebra/EllipseParameters.hh>
#include <numeric/types.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Basic/utility Headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("minimum_bounding_ellipse");

class minimum_bounding_ellipse : public CxxTest::TestSuite {

public:

	void setUp(){

		using namespace numeric;
		using namespace utility;

		core_init();

		// Initialize a vector of 10 random 2D test points
		// (points taken from rand() function in MATLAB)
		test_points_ = vector1< xyzVector< Real > >(10);
		test_points_[1] = xyzVector< Real >(0.8147, 0.9058, 0);
		test_points_[2] = xyzVector< Real >(0.1270, 0.9134, 0);
		test_points_[3] = xyzVector< Real >(0.6324, 0.0975, 0);
		test_points_[4] = xyzVector< Real >(0.2785, 0.5469, 0);
		test_points_[5] = xyzVector< Real >(0.9575, 0.9649, 0);
		test_points_[6] = xyzVector< Real >(0.1576, 0.9706, 0);
		test_points_[7] = xyzVector< Real >(0.9572, 0.4854, 0);
		test_points_[8] = xyzVector< Real >(0.8003, 0.1419, 0);
		test_points_[9] = xyzVector< Real >(0.4218, 0.9157, 0);
		test_points_[10] = xyzVector< Real >(0.7922, 0.9595, 0);

	}

	void tearDown(){

	}

	void test_correct_ellipse_generation() {

		TS_TRACE( "Test correct ellipse generation" );

		using namespace numeric;
		using namespace numeric::linear_algebra;

		// Set tolerance to 0.01
		Real tolerance( 0.01 );

		// Make an ellipse using the 10 random points from MATLAB and a tolerance of 0.01
		EllipseParametersOP ellipse = numeric::linear_algebra::minimum_bounding_ellipse( test_points_, tolerance );

		// Check the accuracy of ellipse parameters
		TS_ASSERT_DELTA( ellipse->center_h(), 0.5797, 0.01 ); // center - h
		TS_ASSERT_DELTA( ellipse->center_k(), 0.6595, 0.01 ); // center - k
		TS_ASSERT_DELTA( ellipse->major_radius(), 0.4737, 0.01 ); // major axis
		TS_ASSERT_DELTA( ellipse->minor_radius(), 0.5608, 0.01 ); // minor axis

		// Check rotation matrix values
		TS_ASSERT_DELTA( ellipse->rotation()(0,0), -0.9667, 0.01 ); // cos(theta)
		TS_ASSERT_DELTA( ellipse->rotation()(0,1), 0.2557, 0.01 ); // -sin(theta)
		TS_ASSERT_DELTA( ellipse->rotation()(1,0), -0.2557, 0.01 ); // sin(theta)
		TS_ASSERT_DELTA( ellipse->rotation()(1,1), -0.9667, 0.01 ); // cos(theta)

	}

	void test_point_in_ellipse() {

		TS_TRACE( "Test point in ellipse" );

		using namespace numeric;
		using namespace numeric::linear_algebra;

		// Set tolerance to 0.01
		Real tolerance( 0.01 );

		// Make an ellipse using the 10 random points from MATLAB and a tolerance of 0.01
		EllipseParametersOP ellipse = numeric::linear_algebra::minimum_bounding_ellipse( test_points_, tolerance );

		// Extrapolate ellipse parameters from above
		Real h( ellipse->center_h() );
		Real k( ellipse->center_k() );
		Real a( ellipse->major_radius() );
		Real b( ellipse->minor_radius() );

		// Going to try applying rotations - will see what happens
		MathMatrix< core::Real > rotation( ellipse->rotation() );

		// Going to try out a few test points based on the toy example
		xyzVector< Real > definitely_out = xyzVector< Real >(100, 100, 0);
		xyzVector< Real > test_point_in = xyzVector< Real >(0.8147, 0.9058, 0);
		xyzVector< Real > test_point_out = xyzVector< Real >(0.1270, 0.9134, 0);
		xyzVector< Real > test_point_on_boundary = xyzVector< Real >(0.4134, 0.5496, 0);

		// hmm I'm wondering if these will pass better if they are rotated? will check ths initially
		TS_ASSERT( !point_in_ellipse( definitely_out, h, k, a, b, rotation ) );
		TS_ASSERT( point_in_ellipse( test_point_in, h, k, a, b, rotation ) );
		TS_ASSERT( !point_in_ellipse( test_point_out, h, k, a, b, rotation ) );
		TS_ASSERT( point_in_ellipse( test_point_on_boundary, h, k, a, b, rotation ) );

	}

private: // data

	utility::vector1< numeric::xyzVector< numeric::Real > > test_points_;

};



