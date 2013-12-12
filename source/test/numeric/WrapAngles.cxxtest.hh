// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  Test suite for numeric/wrap_angles.hh
/// @author Kale Kundert (kale.kundert@ucsf.edu)

// Test headers
#include <cxxtest/TestSuite.h>

// Package Headers
#include <numeric/wrap_angles.hh>
#include <numeric/constants.hh>

using namespace std;
using namespace numeric;
using namespace numeric::constants::r;

class WrapAnglesTest : public CxxTest::TestSuite {

public:

	void test_wrap_360() {
		double inputs[] =  {-450, -360, -270, -180, -90, 0, 90, 180, 270, 360, 450};
		double outputs[] = { 270,    0,   90,  180, 270, 0, 90, 180, 270,   0,  90};

		for (int i=0; i < 9; i++) {
			TS_ASSERT_DELTA(outputs[i], wrap_360(inputs[i]), 1e-10);
		}
	}

	void test_wrap_180() {
		double inputs[] =  {-450, -360, -270, -180, -90, 0, 90,  180, 270, 360, 450};
		double outputs[] = { -90,    0,   90, -180, -90, 0, 90, -180, -90,   0,  90};

		for (int i=0; i < 9; i++) {
			TS_ASSERT_DELTA(outputs[i], wrap_180(inputs[i]), 1e-10);
		}
	}

	void test_wrap_2pi() {
		double inputs[] =  {-2.5*pi, -2.0*pi, -1.5*pi, -1.0*pi, -0.5*pi, 0, 0.5*pi,  1.0*pi,  1.5*pi,  2.0*pi,  2.5*pi};
		double outputs[] = { 1.5*pi,       0,  0.5*pi,  1.0*pi,  1.5*pi, 0, 0.5*pi,  1.0*pi,  1.5*pi,       0,  0.5*pi};

		for (int i=0; i < 9; i++) {
			TS_ASSERT_DELTA(outputs[i], wrap_2pi(inputs[i]), 1e-10);
		}
	}

	void test_wrap_pi() {
		double inputs[] =  {-2.5*pi, -2.0*pi, -1.5*pi, -1.0*pi, -0.5*pi, 0, 0.5*pi,  1.0*pi,  1.5*pi,  2.0*pi,  2.5*pi};
		double outputs[] = {-0.5*pi,       0,  0.5*pi, -1.0*pi, -0.5*pi, 0, 0.5*pi, -1.0*pi, -0.5*pi,       0,  0.5*pi};

		for (int i=0; i < 9; i++) {
			TS_ASSERT_DELTA(outputs[i], wrap_pi(inputs[i]), 1e-10);
		}
	}
};


