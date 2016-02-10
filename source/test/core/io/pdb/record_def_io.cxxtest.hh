// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    test/core/io/pdb/record_def_io.cxxtest.hh
/// @brief   Test suite for database loading of the reference table of PDB record definitions.
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <core/io/pdb/record_def_io.hh>
#include <core/io/pdb/RecordType.hh>
#include <core/io/pdb/Field.hh>

// Basic header
#include <basic/database/open.hh>

// C++ header
#include <utility>
#include <map>


class RecordDefIOTests : public CxxTest::TestSuite {
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
	// Confirm that the reference table of PDB record definitions is loaded correctly.
	void test_read_record_definitions_from_file()
	{
		using namespace core::io::pdb;

		TS_TRACE( "Testing read_record_definitions_from_file() method." );

		std::map< std::string, RecordType > record_type_map;
		record_type_map[ "HEADER" ] = HEADER;
		record_type_map[ "OBSLTE" ] = OBSLTE;
		record_type_map[ "SPLIT" ] = SPLIT;
		record_type_map[ "TITLE" ] = TITLE;
		RecordDef records( read_record_definitions_from_file( "core/io/pdb/fake_pdb_record_defs", record_type_map ) );

		// The 5th record definition in the file should be skipped; it's not in the above map, nor is it a valid key.
		TS_ASSERT_EQUALS( records.size(), 4 );
		TS_ASSERT_EQUALS( records[ HEADER ].size(), 4 );
		TS_ASSERT_EQUALS( records[ OBSLTE ].size(), 4 );
		TS_ASSERT_EQUALS( records[ TITLE ].size(), 3 );
		TS_ASSERT_EQUALS( records[ SPLIT ].size(), 3 );

		TS_ASSERT_EQUALS( records[ HEADER ][ "type" ].value, "" );  // It should be empty when created.
		TS_ASSERT_EQUALS( records[ HEADER ][ "classification" ].start, 11 );
		TS_ASSERT_EQUALS( records[ HEADER ][ "depDate" ].end, 59 );
	}
};  // class RecordDefIOTests
