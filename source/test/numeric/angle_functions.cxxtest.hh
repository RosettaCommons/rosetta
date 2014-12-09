// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/random.cxxtest.hh
/// @brief  test suite for numeric::random
/// @author Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

// Package Headers
#include <numeric/angle.functions.hh>

//#include <vector>
//#include <iostream>
//#include <sstream>

//#include <test/UTracer.hh>

class AngleFunctionTests : public CxxTest::TestSuite
{
public:
	AngleFunctionTests(){}

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	// ------------------------------------------ //
	/// @brief test how Uniform distribution is
	void test_angles() {
		using namespace numeric;
		NumericTraits<double> traits;
		double pi = traits.pi();
		TS_ASSERT_EQUALS(principal_angle(0), 0);
		TS_ASSERT_EQUALS(principal_angle(pi), pi);
		TS_ASSERT_EQUALS(principal_angle(-pi), pi);
		TS_ASSERT_EQUALS(principal_angle(3*pi), -pi);

		TS_ASSERT_EQUALS(principal_angle_radians(0), 0);
		TS_ASSERT_EQUALS(principal_angle_radians(pi), pi);
		TS_ASSERT_EQUALS(principal_angle_radians(-pi), pi);
		TS_ASSERT_EQUALS(principal_angle_radians(3*pi), -pi);

		TS_ASSERT_EQUALS(principal_angle_degrees(0), 0);
		TS_ASSERT_EQUALS(principal_angle_degrees(180), 180);
		TS_ASSERT_EQUALS(principal_angle_degrees(-180), 180);
		TS_ASSERT_EQUALS(principal_angle_degrees(360), 0);

		TS_ASSERT_EQUALS(nonnegative_principal_angle(0), 0);
		TS_ASSERT_EQUALS(nonnegative_principal_angle(pi), pi);
		TS_ASSERT_EQUALS(nonnegative_principal_angle(2*pi), 0);
		TS_ASSERT_EQUALS(nonnegative_principal_angle(3*pi), pi);

		TS_ASSERT_EQUALS(nonnegative_principal_angle_radians(0), 0);
		TS_ASSERT_EQUALS(nonnegative_principal_angle_radians(pi), pi);
		TS_ASSERT_EQUALS(nonnegative_principal_angle_radians(2*pi), 0);
		TS_ASSERT_EQUALS(nonnegative_principal_angle_radians(3*pi), pi);

		TS_ASSERT_EQUALS(nonnegative_principal_angle_degrees(0), 0);
		TS_ASSERT_EQUALS(nonnegative_principal_angle_degrees(180), 180);
		TS_ASSERT_EQUALS(nonnegative_principal_angle_degrees(-180), 180);
		TS_ASSERT_EQUALS(nonnegative_principal_angle_degrees(360), 0);

		TS_ASSERT_EQUALS(nearest_angle(0., pi), 2*pi);
		TS_ASSERT_EQUALS(nearest_angle(1., pi), 1.00);
		TS_ASSERT_EQUALS(nearest_angle(2., pi), 2.00);
		TS_ASSERT_EQUALS(nearest_angle(4, 3), 4.00);

		TS_ASSERT_EQUALS(nearest_angle_radians(pi, 0.), -pi);
		TS_ASSERT_EQUALS(nearest_angle_radians(pi, 1.), pi);
		TS_ASSERT_EQUALS(nearest_angle_radians(pi, 2.), pi);
		TS_ASSERT_EQUALS(nearest_angle(5, 2), 5.00);

		TS_ASSERT_EQUALS(nearest_angle_degrees(0, 180), 0);
		TS_ASSERT_EQUALS(nearest_angle_degrees(90, 180), 90);
		TS_ASSERT_EQUALS(nearest_angle_degrees(180, 180), 180);
		TS_ASSERT_EQUALS(nearest_angle_degrees(90, 30), 90);

	}

};

