// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test/core/chemical/database_io.cxxtest.hh
/// @brief   Test suite for carbohydrate database loading
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/chemical/carbohydrates/database_io.hh>

// Package header
#include <core/chemical/carbohydrates/SugarModificationsNomenclatureTable.hh>

// Project header
#include <core/types.hh>

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

	// Confirm that carbohydrate ring sizes and their 1-letter affixes and morphemes are loaded correctly from the
	// database.
	void test_read_ring_sizes_and_morphemes_from_database_file()
	{
		using namespace std;
		using namespace core::chemical::carbohydrates;

		TS_TRACE( "Testing read_ring_sizes_and_morphemes_from_database_file() method." );

		map< core::Size, pair< char, string > > map( read_ring_sizes_and_morphemes_from_database_file(
			"core/chemical/carbohydrates/ring_size_to_morphemes.map" ) );

		TS_ASSERT_EQUALS( map.size(), 2 );
		TS_ASSERT_EQUALS( map[ 3 ].first, '\0' );  // Make sure 'X' was properly converted to a null char.
		TS_ASSERT_EQUALS( map[ 4 ].second, "ohwow" );

		// Test for bad files.
		TS_TRACE( "An input error should follow:" );
		try {
			read_ring_sizes_and_morphemes_from_database_file(
				"core/chemical/carbohydrates/ring_size_to_morphemes.bad_map" );
		} catch ( utility::excn::EXCN_Base const & e) {
			TS_ASSERT_EQUALS( e.msg().substr( e.msg().find( "ERROR: " ) ),
				"ERROR: read_ring_sizes_and_morphemes_from_database_file: invalid ring size; "
				"rings cannot have less than 3 atoms!\n\n" );
			TS_TRACE( "The above error message was expected." );
		}
	}

	// TODO: Expand this test, once I've finalized what I want/need the table to contain for data.
	// Confirm that the nomenclature data table is loaded correctly from the database.
	void test_read_nomenclature_table_from_database_file()
	{
		using namespace std;
		using namespace core::chemical::carbohydrates;

		TS_TRACE( "Testing read_nomenclature_table_from_database_file() method." );
		SugarModificationsNomenclatureTable table( read_nomenclature_table_from_database_file(
			"core/chemical/carbohydrates/nomenclature.table") );

		TS_ASSERT_EQUALS( table.size(), 8 );  // TEMP
	}
};  // class CarbohydrateDatabaseIOTests
