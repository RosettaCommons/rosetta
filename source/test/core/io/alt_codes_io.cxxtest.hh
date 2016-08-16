// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  test/core/io/pdb/alt_codes_io.cxxtest.hh
/// @brief   Test suite for alternative PDB 3-letter-code database loading
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/io/alt_codes_io.hh>

// Basic header
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// C++ header
#include <utility>

static THREAD_LOCAL basic::Tracer TR("core.io.alt_codes_io.cxxtest");

class AltCodesIOTests : public CxxTest::TestSuite {
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
	// Confirm that alternative PDB 3-letter-codes are loaded correctly from the database.
	void test_read_alternative_3_letter_codes_from_database_file()
	{
		using namespace core::io;

		TR <<  "Testing read_read_alternative_3_letter_codes_from_database_file() method."  << std::endl;

		AltCodeMap alt_codes( read_alternative_3_letter_codes_from_database_file(
			basic::database::full_name( "input_output/3-letter_codes/sentence_case.codes" ) ) );

		TS_ASSERT_EQUALS( alt_codes.size(), 24 );  // 20 canonical and 4 NCAAs
		TS_ASSERT_EQUALS( alt_codes[ "Ala" ].first, "ALA" );
		TS_ASSERT_EQUALS( alt_codes[ "Ala" ].second, "" );
		TS_ASSERT_EQUALS( alt_codes[ "Hcy" ].first, "HCY" );
		TS_ASSERT_EQUALS( alt_codes[ "Hcy" ].second, "HOMOCYSTEINE" );
		TS_ASSERT_EQUALS( alt_codes[ "Val" ].first, "VAL" );
	}
};  // class AltCodesIOTests
