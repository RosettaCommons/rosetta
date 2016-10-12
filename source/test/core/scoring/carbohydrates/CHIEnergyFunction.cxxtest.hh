// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/scoring/CHIEnergyFunction.cxxtest.hh
/// @brief   Test suite for the CHI Energy Functions
/// @author  Labonte <JWLabonte@jhu.edu>


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/scoring/carbohydrates/CHIEnergyFunction.hh>

// Package Header
#include <core/scoring/ScoringManager.hh>

// Project Header
#include <core/types.hh>

// Basic Header
#include <basic/Tracer.hh>


static THREAD_LOCAL basic::Tracer TR( "core.scoring.carbohydrates.CHIEnergyFunction.cxxtest" );


class CHIEnergyFunctionTests : public CxxTest::TestSuite {
public: // Standard methods ///////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		core_init();
	}

	// Destruction
	void tearDown()
	{}


public: // Tests //////////////////////////////////////////////////////////////
	// Confirm that the functions return the correct energies.
	void test_CHI_energy_function()
	{
		using namespace core::scoring::carbohydrates;

		TR <<  "Testing if the CHI Energy Functions return the correct energies."  << std::endl;

		CHIEnergyFunction const & E( core::scoring::ScoringManager::get_instance()->get_CHIEnergyFunction() );
		core::Angle x;

		// All values below calculated with Mathematica.

		x = 0.0;
		TS_ASSERT_DELTA( E( ALPHA_LINKS, x ), 9.42684, 0.002 );
		TS_ASSERT_DELTA( E( BETA_LINKS, x ), 4.99893, 0.002 );
		TS_ASSERT_DELTA( E( _2AX_3EQ_4AX_LINKS, x ), 4.47984, 0.002 );
		TS_ASSERT_DELTA( E( _2EQ_3AX_4EQ_LINKS, x ), 5.49093, 0.002 );
		TS_ASSERT_DELTA( E( ALPHA6_LINKS, x ), 7.4798, 0.002 );
		TS_ASSERT_DELTA( E( BETA6_LINKS, x ), 7.26952, 0.002 );

		x = 60.0;
		TS_ASSERT_DELTA( E( ALPHA_LINKS, x ), 0.431136, 0.002 );
		TS_ASSERT_DELTA( E( BETA_LINKS, x ), 2.85052, 0.002 );
		TS_ASSERT_DELTA( E( _2AX_3EQ_4AX_LINKS, x ), 3.94917, 0.002 );
		TS_ASSERT_DELTA( E( _2EQ_3AX_4EQ_LINKS, x ), 0.980713, 0.002 );
		TS_ASSERT_DELTA( E( ALPHA6_LINKS, x ), 3.61764, 0.002 );
		TS_ASSERT_DELTA( E( BETA6_LINKS, x ), 3.14496, 0.002 );

		x = -60.0;
		TS_ASSERT_DELTA( E( ALPHA_LINKS, x ), 10.3474, 0.002 );
		TS_ASSERT_DELTA( E( BETA_LINKS, x ), 0.119367, 0.002 );
		TS_ASSERT_DELTA( E( _2AX_3EQ_4AX_LINKS, x ), 1.86009, 0.002 );  // x will never be -60 for this function.
		TS_ASSERT_DELTA( E( _2EQ_3AX_4EQ_LINKS, x ), 1.2901, 0.002 );  // x will never be -60 for this function.
		TS_ASSERT_DELTA( E( ALPHA6_LINKS, x ), 67.9606, 0.002 );  // x will never be -60 for this function.
		TS_ASSERT_DELTA( E( BETA6_LINKS, x ), 1.3992, 0.002 );  // x will never be -60 for this function.

		x = -180.0;
		TS_ASSERT_DELTA( E( ALPHA_LINKS, x ), 4.99904, 0.002 );
		TS_ASSERT_DELTA( E( BETA_LINKS, x ), 5.81316, 0.002 );
		TS_ASSERT_DELTA( E( _2AX_3EQ_4AX_LINKS, x ), -0.120704, 0.002 );  // x will never be -180 for this function.
		TS_ASSERT_DELTA( E( _2EQ_3AX_4EQ_LINKS, x ), 1.022, 0.002 );  // x will never be -180 for this function.
		TS_ASSERT_DELTA( E( ALPHA6_LINKS, x ), 0.0000307473, 0.002 );  // x will never be -180 for this function.
		TS_ASSERT_DELTA( E( BETA6_LINKS, x ), 8.07206e-6, 0.002 );  // x will never be -180 for this function.

		x = 180.0;
		TS_ASSERT_DELTA( E( ALPHA_LINKS, x ), 5.00145, 0.002 );
		TS_ASSERT_DELTA( E( BETA_LINKS, x ), 5.81616, 0.002 );
		TS_ASSERT_DELTA( E( _2AX_3EQ_4AX_LINKS, x ), 0.970338, 0.002 );
		TS_ASSERT_DELTA( E( _2EQ_3AX_4EQ_LINKS, x ), 1.12649, 0.002 );
		TS_ASSERT_DELTA( E( ALPHA6_LINKS, x ), 0.017814, 0.002 );
		TS_ASSERT_DELTA( E( BETA6_LINKS, x ), 0.0291781, 0.002 );

		x = 360.0;
		TS_ASSERT_DELTA( E( ALPHA_LINKS, x ), 5.19937e-8, 0.002 );  // x will never be 360 for this function.
		TS_ASSERT_DELTA( E( BETA_LINKS, x ), 14.315, 0.002 );  // x will never be 360 for this function.
		TS_ASSERT_DELTA( E( _2AX_3EQ_4AX_LINKS, x ), 4.47459, 0.002 );
		TS_ASSERT_DELTA( E( _2EQ_3AX_4EQ_LINKS, x ), 5.5771, 0.002 );
		TS_ASSERT_DELTA( E( ALPHA6_LINKS, x ), 7.47945, 0.002 );
		TS_ASSERT_DELTA( E( BETA6_LINKS, x ), 7.22894, 0.002 );
	}

	// Confirm that the function derivatives return the correct values.
	void test_CHI_energy_derivatives()
	{
		using namespace core::scoring::carbohydrates;

		TR <<  "Testing if the CHI Energy Function derivatives return the correct values."  << std::endl;

		CHIEnergyFunction const & E( core::scoring::ScoringManager::get_instance()->get_CHIEnergyFunction() );
		core::Angle x;

		// All values below calculated with Mathematica.

		x = 0.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ALPHA_LINKS, x ), -0.117671, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( BETA_LINKS, x ), 0.0251501, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( _2AX_3EQ_4AX_LINKS, x ), 0.00994998, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( _2EQ_3AX_4EQ_LINKS, x ), 0.0000290124, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( ALPHA6_LINKS, x ), -0.09549, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( BETA6_LINKS, x ), 0.025375, 0.002 );

		x = 60.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ALPHA_LINKS, x ), -0.0589094, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( BETA_LINKS, x ), 0.0353785, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( _2AX_3EQ_4AX_LINKS, x ), 0.0106208, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( _2EQ_3AX_4EQ_LINKS, x ), -0.0466017, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( ALPHA6_LINKS, x ), -0.0895345, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( BETA6_LINKS, x ), -0.059447, 0.002 );

		x = -60.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ALPHA_LINKS, x ), -0.00620749, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( BETA_LINKS, x ), 0.0305844, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( _2AX_3EQ_4AX_LINKS, x ), 0.0516056, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( _2EQ_3AX_4EQ_LINKS, x ), 0.0251424, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( ALPHA6_LINKS, x ), 0.0678384, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( BETA6_LINKS, x ), 0.0723912, 0.002 );

		x = -180.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ALPHA_LINKS, x ), 0.00656956, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( BETA_LINKS, x ), -0.136412, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( _2AX_3EQ_4AX_LINKS, x ), 0.000365645, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( _2EQ_3AX_4EQ_LINKS, x ), 1.36255e-11, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( ALPHA6_LINKS, x ), 7.45748e-6, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( BETA6_LINKS, x ), 1.20531e-6, 0.002 );

		x = 180.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ALPHA_LINKS, x ), 0.0741118, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( BETA_LINKS, x ), -0.0133324, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( _2AX_3EQ_4AX_LINKS, x ), -0.043664, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( _2EQ_3AX_4EQ_LINKS, x ), 0.0445609, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( ALPHA6_LINKS, x ), -0.0134757, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( BETA6_LINKS, x ), -0.0417507, 0.002 );

		x = 360.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ALPHA_LINKS, x ), -1.1386e-8, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( BETA_LINKS, x ), -0.217421, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( _2AX_3EQ_4AX_LINKS, x ), 0.0109582, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( _2EQ_3AX_4EQ_LINKS, x ), -0.0115829, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( ALPHA6_LINKS, x ), -0.0636112, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( BETA6_LINKS, x ), 0.0381025, 0.002 );
	}

	void test_CHI_sampling_data()
	{
		using core::scoring::ScoringManager;
		using namespace core::scoring::carbohydrates;

		CHIEnergyFunction const & sugar_bb =
			ScoringManager::get_instance()->get_CHIEnergyFunction( true /*setup for scoring*/, .1 );


		TS_ASSERT( sugar_bb.sampling_data_setup() );
		for ( core::Size i = 1; i <= 4; ++i ) {
			TR << "Linkage type: " << i << std::endl;


			CHIEnergyFunctionLinkageType link_type = static_cast< CHIEnergyFunctionLinkageType >( i );

			TS_ASSERT( sugar_bb.sampling_data_setup( link_type ) );
			TS_ASSERT_THROWS_NOTHING( sugar_bb.get_chi_sampling_data( link_type )  );
			CHIDihedralSamplingData const & sampling_data = sugar_bb.get_chi_sampling_data( link_type );

			TR << "PROB Size" << sampling_data.probabilities.size() << std::endl;
		}
	}
};
