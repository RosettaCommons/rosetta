// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/random/reservoir_sample.cxxtest.hh
/// @brief
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

// Package Headers
#include <numeric/random/random.hh>
#include <numeric/random/reservoir_sample.hh>

#include <vector>
#include <iostream>
#include <sstream>

#include <test/UTracer.hh>


class ReservoirSamplingTests : public CxxTest::TestSuite {
public:
	ReservoirSamplingTests() {
		numeric::random::rg().set_seed(	"mt19937", 999 );
	}

	// Shared initialization goes here.
	void setUp() {}

	// Shared finalization goes here.
	void tearDown() {}

	void test_reservoir_sample() {
		using numeric::Size;
		using numeric::Real;
		using utility::vector1;
		using numeric::random::reservoir_sample;

		vector1< Size > original_vals;
		for ( Size ii = 1; ii <= 10; ++ii ) {
			original_vals.push_back( ii );
		}

		vector1< Size > counts( original_vals.size(), 0 );
		Size const n_wanted( 2 );
		Size const total( 10000 );
		for ( Size jj = 1; jj <= total; ++jj ) {
			vector1< Size > sample = reservoir_sample< Size >(
				original_vals, n_wanted, numeric::random::rg()
			);

			for ( Size ii = 1; ii <= sample.size(); ++ii ) {
				counts[ sample[ ii ] ]++;
			}

			TS_ASSERT_EQUALS( sample.size(), n_wanted );
		}

		Real const expected(
			static_cast< Real >(n_wanted) / static_cast< Real >(original_vals.size())
		);
		for ( Size jj = 1; jj <= counts.size(); ++jj ) {
			Real const observed(
				static_cast< Real >(counts[jj]) / static_cast< Real >(total)
			);
			TS_ASSERT_DELTA( observed, expected, 1e-2 );
		}

		//using numeric::Real;
		//for ( Size ii = 1; ii <= counts.size(); ++ii ) {
		//	Real const percent(
		//		static_cast< Real >(counts[ii]) / static_cast< Real >(total)
		//	);
		//	std::cout << ii << ' ' << counts[ii] << ' ' << percent << std::endl;
		//}

	}

};
