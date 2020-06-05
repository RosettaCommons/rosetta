// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/random/DistributionSampler.cxxtest.hh
/// @brief test suite for numeric/random/DistributionSampler
/// @author Christopher Miles (cmiles@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

// C/C++ headers
#include <cmath>
#include <iostream>

// External headers

// Project headers
#include <core/types.hh>
#include <numeric/random/DistributionSampler.hh>

#define TOLERANCE 0.1
#define NUM_SAMPLES 1000000

namespace {

using boost::math::normal;
using numeric::random::DistributionSampler;

class DistributionSamplerTest : public CxxTest::TestSuite {
public:
	// If we draw a sufficient number of samples, can we recover the underlying
	// distribution's parameters?
	void test_normal() {

		// define a normal distribution with mean 0 and std 1
		normal dist(0, 1);
		DistributionSampler<normal> sampler(dist);

		// draw samples from the distribution
		// This used to use boost accumulators, but ICC seems to have issues with them.
		utility::vector1< core::Real > samples;
		samples.reserve( NUM_SAMPLES );
		for ( int i = 0; i < NUM_SAMPLES; ++i ) {
			samples.push_back( sampler.sample() );
		}

		// calculate population mean and std
		core::Real population_mean = std::accumulate( samples.begin(), samples.end(), 0 ) / NUM_SAMPLES;

		core::Real population_std  = 0;
		for ( core::Real x: samples ) {
			core::Real const d = x - population_mean;
			population_std += d * d;
		}
		population_std = std::sqrt( population_std / NUM_SAMPLES );

		TS_ASSERT_DELTA(dist.mean(), population_mean, TOLERANCE);
		TS_ASSERT_DELTA(dist.standard_deviation(), population_std, TOLERANCE);
	}
};

}  // anonymous namespace
