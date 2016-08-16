// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/color_util.cxxtest.hh
/// @author Sam DeLuca
/// @brief test suite for colorspace utility functions

// Test colors
#include <cxxtest/TestSuite.h>
#include <numeric/color_util.hh>
#include <numeric/xyzVector.hh>

class ColorUtilTests : public CxxTest::TestSuite {

public:

	platform::Real delta;
	void setUp()
	{
		delta = 0.01;
	}

	void tearDown()
	{

	}

	void test_rgb_hsv()
	{
		numeric::xyzVector<platform::Real> test_white_hsv(numeric::rgb_to_hsv(1.0,1.0,1.0));
		numeric::xyzVector<platform::Real> real_white_hsv(0,0,1.00);
		TS_ASSERT(test_white_hsv == real_white_hsv);


		numeric::xyzVector<platform::Real> test_black_hsv(numeric::rgb_to_hsv(0,0,0));
		numeric::xyzVector<platform::Real> real_black_hsv(0,0,0);
		TS_ASSERT(test_black_hsv == real_black_hsv);

		numeric::xyzVector<platform::Real> test_grey_hsv(numeric::rgb_to_hsv(0.39,0.39,0.39));
		numeric::xyzVector<platform::Real> real_grey_hsv(0,0,0.39);
		TS_ASSERT(test_grey_hsv == real_grey_hsv);

		numeric::xyzVector<platform::Real> test_red_hsv(numeric::rgb_to_hsv(1.0,0,0));
		numeric::xyzVector<platform::Real> real_red_hsv(0,1.0,1.0);
		TS_ASSERT(test_red_hsv == real_red_hsv);

		numeric::xyzVector<platform::Real> test_green_hsv(numeric::rgb_to_hsv(0,1.0,0));
		numeric::xyzVector<platform::Real> real_green_hsv(120,1.0,1.0);
		TS_ASSERT(test_green_hsv == real_green_hsv);

		numeric::xyzVector<platform::Real> test_blue_hsv(numeric::rgb_to_hsv(0,0,1.0));
		numeric::xyzVector<platform::Real> real_blue_hsv(240,1.0,1.0);
		TS_ASSERT(test_blue_hsv == real_blue_hsv);

		numeric::xyzVector<platform::Real> test_random_hsv(numeric::rgb_to_hsv(0.14,0.25,0.58));
		numeric::xyzVector<platform::Real> real_random_hsv(225,0.76,0.58);

		TS_ASSERT_DELTA(test_random_hsv.x(),real_random_hsv.x(),delta);
		TS_ASSERT_DELTA(test_random_hsv.y(),real_random_hsv.y(),delta);
		TS_ASSERT_DELTA(test_random_hsv.z(),real_random_hsv.z(),delta);

	}

	void test_hsv_rgb()
	{
		numeric::xyzVector<platform::Real> test_white_rgb(numeric::hsv_to_rgb(0,0,1.00));
		numeric::xyzVector<platform::Real> real_white_rgb(1.0,1.0,1.0);
		TS_ASSERT(test_white_rgb == real_white_rgb);

		numeric::xyzVector<platform::Real> test_black_rgb(numeric::hsv_to_rgb(0,0,0));
		numeric::xyzVector<platform::Real> real_black_rgb(0,0,0);
		TS_ASSERT(test_black_rgb == real_black_rgb);

		numeric::xyzVector<platform::Real> test_grey_rgb(numeric::hsv_to_rgb(0,0,0.39));
		numeric::xyzVector<platform::Real> real_grey_rgb(0.39,0.39,0.39);
		TS_ASSERT(test_grey_rgb == real_grey_rgb);

		numeric::xyzVector<platform::Real> test_red_rgb(numeric::hsv_to_rgb(0,1.0,1.0));
		numeric::xyzVector<platform::Real> real_red_rgb(1.0,0,0);
		TS_ASSERT(test_red_rgb == real_red_rgb);

		numeric::xyzVector<platform::Real> test_green_rgb(numeric::hsv_to_rgb(120,1.0,1.0));
		numeric::xyzVector<platform::Real> real_green_rgb(0.0,1.0,0.0);
		TS_ASSERT(test_green_rgb == real_green_rgb);

		numeric::xyzVector<platform::Real> test_blue_rgb(numeric::hsv_to_rgb(240,1.0,1.0));
		numeric::xyzVector<platform::Real> real_blue_rgb(0,0,1.0);
		TS_ASSERT(test_blue_rgb == real_blue_rgb);

		numeric::xyzVector<platform::Real> test_random_rgb(numeric::hsv_to_rgb(255,0.77,0.58));
		numeric::xyzVector<platform::Real> real_random_rgb(0.14,0.25,0.58);


	}
};
