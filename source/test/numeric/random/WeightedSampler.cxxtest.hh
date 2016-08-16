// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/random/WeightedSampler.cxxtest.hh
/// @brief  test suite for numeric::random::WeightedSampler
/// @author Colin A. Smith

// Test headers
#include <cxxtest/TestSuite.h>

// Package Headers
#include <numeric/random/WeightedSampler.hh>

#include <numeric/constants.hh>

#include <vector>
#include <iostream>
#include <sstream>


class WeightedSamplerTests : public CxxTest::TestSuite
{
public:

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// ------------------------------------------ //
	/// @brief test WeightedSampler
	void test_WeightedSampler() {

		numeric::random::WeightedSampler ws1;

		// simple case, two elements with identical weights

		ws1.clear();
		ws1.add_weight(1);
		ws1.add_weight(1);

		TS_ASSERT_EQUALS(ws1.size(), (unsigned)2);
		TS_ASSERT_EQUALS(ws1.random_sample(0), (unsigned)1);
		TS_ASSERT_EQUALS(ws1.random_sample(0), (unsigned)1);
		TS_ASSERT_EQUALS(ws1.random_sample(0.25), (unsigned)1);
		TS_ASSERT_EQUALS(ws1.random_sample(0.5), (unsigned)1);
		TS_ASSERT_EQUALS(ws1.random_sample(0.75), (unsigned)2);
		TS_ASSERT_EQUALS(ws1.random_sample(1), (unsigned)2);

		// more complex case, 0 weighted elements surrounding the identical weights

		ws1.clear();
		ws1.add_weight(0);
		ws1.add_weight(1);
		ws1.add_weight(0);
		ws1.add_weight(1);
		ws1.add_weight(0);

		TS_ASSERT_EQUALS(ws1.size(), (unsigned)5);
		TS_ASSERT_EQUALS(ws1.random_sample(0), (unsigned)2);
		TS_ASSERT_EQUALS(ws1.random_sample(0.25), (unsigned)2);
		TS_ASSERT_EQUALS(ws1.random_sample(0.5), (unsigned)2);
		TS_ASSERT_EQUALS(ws1.random_sample(0.75), (unsigned)4);
		TS_ASSERT_EQUALS(ws1.random_sample(1), (unsigned)4);

		// test size constructor

		numeric::random::WeightedSampler ws2(5);
		ws2.set_weight(2, 1);
		ws2.set_weight(4, 1);

		TS_ASSERT_EQUALS(ws2.size(), (unsigned)5);
		TS_ASSERT_EQUALS(ws2.random_sample(0), (unsigned)2);
		TS_ASSERT_EQUALS(ws2.random_sample(0.25), (unsigned)2);
		TS_ASSERT_EQUALS(ws2.random_sample(0.5), (unsigned)2);
		TS_ASSERT_EQUALS(ws2.random_sample(0.75), (unsigned)4);
		TS_ASSERT_EQUALS(ws2.random_sample(1), (unsigned)4);

		// test vector constructor

		utility::vector1<numeric::Real> weights;
		weights.push_back(1);
		weights.push_back(1);
		numeric::random::WeightedSampler ws3(weights);

		TS_ASSERT_EQUALS(ws3.size(), (unsigned)2);
		TS_ASSERT_EQUALS(ws3.random_sample(0), (unsigned)1);
		TS_ASSERT_EQUALS(ws3.random_sample(0), (unsigned)1);
		TS_ASSERT_EQUALS(ws3.random_sample(0.25), (unsigned)1);
		TS_ASSERT_EQUALS(ws3.random_sample(0.5), (unsigned)1);
		TS_ASSERT_EQUALS(ws3.random_sample(0.75), (unsigned)2);
		TS_ASSERT_EQUALS(ws3.random_sample(1), (unsigned)2);
	}
};

