// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/interpolate.cxxtest.hh
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <cxxtest/TestSuite.h>

// C/C++ headers
#include <iterator>
#include <vector>

// Utility headers
#include <utility/vector1.hh>

// Project headers
#include <numeric/interpolate.hh>

#define MARGIN 0.01

using namespace std;

class InterpolateTest : public CxxTest::TestSuite {
 public:
  void test_linear_ascending() {
    double a = 0;
    double b = 1;
    unsigned num_stages = 2;
    TS_ASSERT_DELTA(0, numeric::linear_interpolate(a, b, 0, num_stages), MARGIN);
    TS_ASSERT_DELTA(0.5, numeric::linear_interpolate(a, b, 1, num_stages), MARGIN);
    TS_ASSERT_DELTA(1, numeric::linear_interpolate(a, b, 2, num_stages), MARGIN);
  }

  void test_linear_descending() {
    double a = 1;
    double b = 0;
    unsigned num_stages = 2;
    TS_ASSERT_DELTA(1, numeric::linear_interpolate(a, b, 0, num_stages), MARGIN);
    TS_ASSERT_DELTA(0.5, numeric::linear_interpolate(a, b, 1, num_stages), MARGIN);
    TS_ASSERT_DELTA(0, numeric::linear_interpolate(a, b, 2, num_stages), MARGIN);
  }
};
