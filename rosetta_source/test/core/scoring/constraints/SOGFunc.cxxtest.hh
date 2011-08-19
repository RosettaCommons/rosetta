// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/SOGFunc.cxxtest.hh
/// @brief  test suite for SOGFunc function
/// @author James Thompson

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <core/scoring/constraints/SOGFunc.hh>
#include <core/scoring/constraints/SOGFunc.fwd.hh>

#include <core/types.hh>


class SOGFuncTests : public CxxTest::TestSuite {

public:
	SOGFuncTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {}

	void test_func() {
		using core::Size;
		using core::Real;
		using utility::vector1;
		using namespace core::scoring::constraints;

		vector1< Real > means, sdevs, weights;
		means.push_back( 0.5 );
		sdevs.push_back( 1.0 );
		weights.push_back( 0.3333 );

		means.push_back( 1.5 );
		sdevs.push_back( 0.5 );
		weights.push_back( 0.3333 );

		means.push_back( 5.0 );
		sdevs.push_back( 2.0 );
		weights.push_back( 0.3333 );

		SOGFuncOP func( new SOGFunc( means, sdevs, weights ) );

    float const TOLERATED_ERROR = 0.001;
		core::Real const start(  2.0 );
		core::Real const end  ( 20.0 );
		core::Real const res  (  0.5 );

		core::Real func_values[] = {
			1.487, 2.4719, 3.01347, 2.96131, 2.83084, 2.74136, 2.71072, 2.74204,
			2.8358, 2.99205, 3.2108, 3.49205, 3.8358, 4.24205, 4.7108, 5.24205,
			5.8358, 6.49205, 7.2108, 7.99205, 8.8358, 9.74205, 10.7108, 11.742,
			12.8358, 13.992, 15.2108, 16.492, 17.8358, 19.242, 20.7108, 22.242,
			23.8358, 25.492, 27.2108, 28.992
		};

		core::Real dfunc_values[] = {
			1.642, 1.906, 0.248, -0.264, -0.231, -0.122, 0.000, 0.125, 0.250, 0.375,
			0.500, 0.625, 0.750, 0.875, 1.000, 1.125, 1.250, 1.375, 1.500, 1.625, 1.750,
			1.875, 2.000, 2.125, 2.250, 2.375, 2.500, 2.625, 2.750, 2.875, 3.000, 3.125,
			3.250, 3.375, 3.500, 3.625, 3.750
		};

		core::Size const nsteps( core::Size( ( end - start ) / res ) );
		for ( core::Size i = 0; i < nsteps; ++i ) {
			core::Real r = start + (i * res);
			TS_ASSERT_DELTA( func->func(r),   func_values[i], TOLERATED_ERROR );
			TS_ASSERT_DELTA( func->dfunc(r), dfunc_values[i], TOLERATED_ERROR );
			TS_ASSERT_DELTA( func->dfunc(r), func->estimate_dfunc(r), TOLERATED_ERROR );
		}
	} // test_func
};
