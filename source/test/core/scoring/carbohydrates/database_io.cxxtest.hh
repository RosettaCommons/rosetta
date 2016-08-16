// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  test/core/scoring/database_io.cxxtest.hh
/// @brief   Test suite for carbohydrate scoring database loading
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/scoring/carbohydrates/database_io.hh>

// Project header
#include <core/types.hh>

// Utility header
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ header
#include <map>

static THREAD_LOCAL basic::Tracer TR("core.scoring.carbohydrates.database_io.cxxtest");

class CarbohydrateScoringDatabaseIOTests : public CxxTest::TestSuite {
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
	// Confirm that Gaussian parameters are loaded correctly from the database.
	void test_read_Gaussian_parameters_from_database_file()
	{
		using namespace std;
		using namespace utility;
		using namespace core::scoring::carbohydrates;

		TR <<  "Testing read_Gaussian_parameters_from_database_file() method."  << std::endl;

		map< char, vector1< core::Real > > parameters(
			read_Gaussian_parameters_from_database_file( "core/scoring/carbohydrates/dummy_params.data" ) );

		TS_ASSERT_EQUALS( parameters.size(), 6);
		TS_ASSERT_EQUALS( parameters[ 'f' ].size(), 1 );
		TS_ASSERT_EQUALS( parameters[ 'o' ].size(), 2 );
		TS_ASSERT_EQUALS( parameters[ 'O' ].size(), 3 );
		TS_ASSERT_EQUALS( parameters[ 'b' ].size(), 4 );
		TS_ASSERT_EQUALS( parameters[ 'a' ].size(), 5 );
		TS_ASSERT_EQUALS( parameters[ 'r' ].size(), 6 );

		TS_ASSERT_EQUALS( parameters[ 'r' ][ 1 ], -6.78e-1 );
		TS_ASSERT_EQUALS( parameters[ 'r' ][ 2 ], 9.01 );
		TS_ASSERT_EQUALS( parameters[ 'r' ][ 3 ], -2.34 );
		TS_ASSERT_EQUALS( parameters[ 'r' ][ 4 ], 5.67 );
		TS_ASSERT_EQUALS( parameters[ 'r' ][ 5 ], -8.90 );
		TS_ASSERT_EQUALS( parameters[ 'r' ][ 6 ], 1.23 );
	}
};  // class CarbohydrateScoringDatabaseIOTests
