// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	 test/core/chemical/database_io.cxxtest.hh
/// @brief   Test suite for carbohydrate database loading
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/chemical/carbohydrates/database_io.hh>

// C++ header
#include <map>


class CarbohydrateDatabaseIOTests : public CxxTest::TestSuite {
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
	// Confirm that carbohydrate 3-letter codes and roots are loaded correctly from the database.
	void test_read_codes_and_roots_from_database_file()
	{
		using namespace std;
		using namespace core::chemical::carbohydrates;

		TS_TRACE( "Testing read_codes_and_roots_from_database_file() method." );

		map< string, string > map(
				read_codes_and_roots_from_database_file( "core/chemical/carbohydrates/codes_to_roots.map" ) );

		TS_ASSERT_EQUALS( map.size(), 3 );
	}
};  // class CarbohydrateDatabaseIOTests
