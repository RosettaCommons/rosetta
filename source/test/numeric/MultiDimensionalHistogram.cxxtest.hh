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
#include <numeric/MultiDimensionalHistogram.hh>

#include <numeric/constants.hh>

#include <vector>
#include <iostream>
#include <sstream>


class MultiDimensionalHistogramTests : public CxxTest::TestSuite
{
public:

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// ------------------------------------------ //
	/// @brief test MultiDimensionalHistogram
	void test_MultiDimensionalHistogram() {

		numeric::MultiDimensionalHistogram mdhist(1);

		// standard 1D histogram
		mdhist.label("1D");
		mdhist.set_dimension(1, 2, 0, numeric::constants::r::pi_2, "D1");

		// underflow
		mdhist.record(-0.0000001);
		TS_ASSERT_EQUALS((signed)mdhist.counts()[0], 1);

		// start
		mdhist.record(0);
		TS_ASSERT_EQUALS((signed)mdhist.counts()[1], 1);

		// bin boundary rounds up
		mdhist.record(numeric::constants::r::pi);
		TS_ASSERT_EQUALS((signed)mdhist.counts()[2], 1);

		// end
		mdhist.record(numeric::constants::r::pi_2);
		TS_ASSERT_EQUALS((signed)mdhist.counts()[2], 2);

		// overflow
		mdhist.record(numeric::constants::r::pi_2+0.0000001);
		TS_ASSERT_EQUALS((signed)mdhist.counts()[3], 1);

		// middle
		mdhist.record(0.5);
		TS_ASSERT_EQUALS((signed)mdhist.counts()[1], 2);


		// default 1D histogram
		numeric::MultiDimensionalHistogram mdhist2;

		// underflow
		mdhist2.record(-1);
		TS_ASSERT_EQUALS((signed)mdhist2.counts()[0], 1);

		// zero
		mdhist2.record(0);
		TS_ASSERT_EQUALS((signed)mdhist2.counts()[1], 1);

		// overflow
		mdhist2.record(1);
		TS_ASSERT_EQUALS((signed)mdhist2.counts()[2], 1);


		// 2D histogram
		numeric::MultiDimensionalHistogram mdhist3(2);
		mdhist3.label("2D");
		mdhist3.set_dimension(1, 2, 0, numeric::constants::r::pi_2, "D1");
		mdhist3.set_dimension(2, 2, 0, numeric::constants::r::pi_2, "D2");
		utility::vector1<numeric::Real> values;
		values.resize(2);

		// 0,0
		mdhist3.record(values);
		TS_ASSERT_EQUALS((signed)mdhist3.counts()[1+1*4], 1);

		// 2*pi,2*pi
		values[1] = numeric::constants::r::pi_2;
		values[2] = numeric::constants::r::pi_2;
		mdhist3.record(values);
		TS_ASSERT_EQUALS((signed)mdhist3.counts()[2+2*4], 1);

		// underflow,overflow
		values[1] = -0.001;
		values[2] = numeric::constants::r::pi_2+0.001;
		mdhist3.record(values);
		TS_ASSERT_EQUALS((signed)mdhist3.counts()[0+3*4], 1);

		//std::cout << mdhist << std::endl << mdhist2 << std::endl << mdhist3;
	}
};

