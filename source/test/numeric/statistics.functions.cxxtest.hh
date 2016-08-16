// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/statistics.functions.cxxtest.hh
/// @brief  test suite for statistics_functions
/// @author James Thompson

#include <cxxtest/TestSuite.h>

#include <utility/vector1.hh>
#include <numeric/types.hh>
#include <numeric/statistics/functions.hh>

class StatisticsFunctionsTests : public CxxTest::TestSuite {

	public:

	double delta_percent;

	void setUp() {
		delta_percent = 0.0001;
	}

	void tearDown() {}

	void test_kl_divergence() {
		using numeric::Real;
		using utility::vector1;
		using numeric::statistics::kl_divergence;

		vector1< Real > prior, posterior;
	  prior.push_back( 0.05467002 );
		prior.push_back( 0.08065691 );
		prior.push_back( 0.10648267 );
		prior.push_back( 0.12579441 );
		prior.push_back( 0.13298076 );
		prior.push_back( 0.12579441 );
		prior.push_back( 0.10648267 );
		prior.push_back( 0.08065691 );
		prior.push_back( 0.05467002 );
		prior.push_back( 0.03315905 );
		posterior.push_back( 1.338302e-04 );
		posterior.push_back( 4.431848e-03 );
		posterior.push_back( 5.399097e-02 );
		posterior.push_back( 2.419707e-01 );
		posterior.push_back( 3.989423e-01 );
		posterior.push_back( 2.419707e-01 );
		posterior.push_back( 5.399097e-02 );
		posterior.push_back( 4.431848e-03 );
		posterior.push_back( 1.338302e-04 );
		posterior.push_back( 1.486720e-06 );
		TS_ASSERT_DELTA( kl_divergence( prior, posterior ), 1.2914, delta_percent );
	}
};
