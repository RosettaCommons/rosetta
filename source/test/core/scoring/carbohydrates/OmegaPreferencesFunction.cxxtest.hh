// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/scoring/OmegaPreferencesFunction.cxxtest.hh
/// @brief   Test suite for the omega-preferences functions
/// @author  Labonte <JWLabonte@jhu.edu>


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/scoring/carbohydrates/OmegaPreferencesFunction.hh>

// Package Header
#include <core/scoring/ScoringManager.hh>

// Project Header
#include <core/types.hh>

// Basic Header
#include <basic/Tracer.hh>


static basic::Tracer TR( "core.scoring.carbohydrates.OmegaPreferencesFunction.cxxtest" );


class OmegaPreferencesFunctionTests : public CxxTest::TestSuite {
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
	void test_omega_preferences_function()
	{
		using namespace core::scoring::carbohydrates;

		TR <<  "Testing whether the omega-preferences functions return the correct energies."  << std::endl;

#ifdef MULTI_THREADED
		try {
			core::scoring::ScoringManager::get_instance()->get_OmegaPreferencesFunction();
		} catch (utility::excn::Exception& excn )  {
			TR << excn.msg() << std::endl;
			std::string expected( "ERROR: Error in ScoringManager: the carbohydrate OmegaPreferencesFunction is fundamentally not threadsafe, and cannot be used in a multithreaded environment.  Please contact Jason Labonte (JWLabonte@jhu.edu) to complain about this." );
			TS_ASSERT_EQUALS( excn.msg().substr( excn.msg().find( "ERROR: " ), expected.size() ), expected );
		}
		return;
#endif

		OmegaPreferencesFunction const & E( core::scoring::ScoringManager::get_instance()->get_OmegaPreferencesFunction() );
		core::Angle x;

		// All values below calculated with Mathematica or mentally.

		x = -60.0;
		TS_ASSERT_DELTA( E( ANTI, x ), 36.0, 0.002 );  // x will never be -60 for this function.
		TS_ASSERT_DELTA( E( GAUCHE_EFFECT, x ), 36.21, 0.002 );  // x will never be -60 for this function.

		x = 0.0;
		TS_ASSERT_DELTA( E( ANTI, x ), 9.0, 0.002 );
		TS_ASSERT_DELTA( E( GAUCHE_EFFECT, x ), 9.21, 0.002 );

		x = 30.0;
		TS_ASSERT_DELTA( E( ANTI, x ), 2.25, 0.002 );
		TS_ASSERT_DELTA( E( GAUCHE_EFFECT, x ), 2.46, 0.002 );

		x = 60.0;
		TS_ASSERT_DELTA( E( ANTI, x ), 0.0, 0.002 );
		TS_ASSERT_DELTA( E( GAUCHE_EFFECT, x ), 0.21, 0.002 );

		x = 120.0;
		TS_ASSERT_DELTA( E( ANTI, x ), 9.0, 0.002 );
		TS_ASSERT_DELTA( E( GAUCHE_EFFECT, x ), 9.21, 0.002 );

		x = 180.0;
		TS_ASSERT_DELTA( E( ANTI, x ), 0.3, 0.002 );
		TS_ASSERT_DELTA( E( GAUCHE_EFFECT, x ), 1.39, 0.002 );

		x = 240.0;
		TS_ASSERT_DELTA( E( ANTI, x ), 9.3, 0.002 );
		TS_ASSERT_DELTA( E( GAUCHE_EFFECT, x ), 10.39, 0.002 );

		x = 300.0;
		TS_ASSERT_DELTA( E( ANTI, x ), 1.0, 0.002 );
		TS_ASSERT_DELTA( E( GAUCHE_EFFECT, x ), 0.0, 0.002 );

		x = 360.0;
		TS_ASSERT_DELTA( E( ANTI, x ), 10.0, 0.002 );
		TS_ASSERT_DELTA( E( GAUCHE_EFFECT, x ), 9.0, 0.002 );
	}

	// Confirm that the function derivatives return the correct values.
	void test_omega_preferences_derivatives()
	{
		using namespace core::scoring::carbohydrates;

		TR <<  "Testing whether the omega-preferences function derivatives return the correct values."  << std::endl;

#ifdef MULTI_THREADED
		try {
			core::scoring::ScoringManager::get_instance()->get_OmegaPreferencesFunction();
		} catch (utility::excn::Exception& excn )  {
			TR << excn.msg() << std::endl;
			std::string expected( "ERROR: Error in ScoringManager: the carbohydrate OmegaPreferencesFunction is fundamentally not threadsafe, and cannot be used in a multithreaded environment.  Please contact Jason Labonte (JWLabonte@jhu.edu) to complain about this." );
			TS_ASSERT_EQUALS( excn.msg().substr( excn.msg().find( "ERROR: " ), expected.size() ), expected );
		}
		return;
#endif

		OmegaPreferencesFunction const & E( core::scoring::ScoringManager::get_instance()->get_OmegaPreferencesFunction() );
		core::Angle x;

		// All values below calculated with Mathematica or mentally.ygeyl01

		x = -60.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ANTI, x ), -0.6, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( GAUCHE_EFFECT, x ), -0.6, 0.002 );

		x = 0.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ANTI, x ), -0.3, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( GAUCHE_EFFECT, x ), -0.3, 0.002 );

		x = 30.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ANTI, x ), -0.15, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( GAUCHE_EFFECT, x ), -0.15, 0.002 );

		x = 60.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ANTI, x ), 0.0, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( GAUCHE_EFFECT, x ), 0.0, 0.002 );

		x = 120.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ANTI, x ), 0.3, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( GAUCHE_EFFECT, x ), 0.3, 0.002 );

		x = 180.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ANTI, x ), 0.0, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( GAUCHE_EFFECT, x ), 0.0, 0.002 );

		x = 240.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ANTI, x ), 0.3, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( GAUCHE_EFFECT, x ), 0.3, 0.002 );

		x = 300.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ANTI, x ), 0.0, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( GAUCHE_EFFECT, x ), 0.0, 0.002 );

		x = 360.0;
		TS_ASSERT_DELTA( E.evaluate_derivative( ANTI, x ), 0.3, 0.002 );
		TS_ASSERT_DELTA( E.evaluate_derivative( GAUCHE_EFFECT, x ), 0.3, 0.002 );
	}
};
