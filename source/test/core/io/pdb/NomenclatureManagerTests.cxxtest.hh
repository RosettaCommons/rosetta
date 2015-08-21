// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  test/core/io/pdb/NomenclatureManagerTests.cxxtest.hh
/// @brief   Test suite for the NomenclatureManager singleton
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/io/pdb/NomenclatureManager.hh>

// Utility header
#include <utility/vector1.hh>

// Basic headers
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// C++ header
#include <string>


class NomenclatureManagerTests : public CxxTest::TestSuite {
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
	// Confirm that the NomenclatureManager hands back the expected name3 and ResidueType base names.
	void test_rosetta_names_from_pdb_code()
	{
		using namespace std;
		using namespace basic::options;
		using namespace core::io::pdb;

		TS_TRACE( "Testing rosetta_names_from_pdb_code() static method with alternative 3-letter codes not provided." );

		pair< string, string > residue_empty(
			NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "" ) );
		pair< string, string > residue_sentence_case(
			NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "Ala" ) );
		pair< string, string > residue_all_caps(
			NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "ALA" ) );
		pair< string, string > residue_special(
			NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "Hcy" ) );
		pair< string, string > residue_bogus(
			NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "WOO!" ) );

		TS_ASSERT_EQUALS( residue_empty.first, "" );
		TS_ASSERT_EQUALS( residue_empty.second, "" );
		TS_ASSERT_EQUALS( residue_sentence_case.first, "Ala" );
		TS_ASSERT_EQUALS( residue_sentence_case.second, "" );
		TS_ASSERT_EQUALS( residue_all_caps.first, "ALA" );
		TS_ASSERT_EQUALS( residue_all_caps.second, "" );
		TS_ASSERT_EQUALS( residue_special.first, "Hcy" );
		TS_ASSERT_EQUALS( residue_special.second, "" );
		TS_ASSERT_EQUALS( residue_bogus.first, "WOO!" );
		TS_ASSERT_EQUALS( residue_bogus.second, "" );

		TS_TRACE( "Testing rosetta_names_from_pdb_code() static method with alternative 3-letter codes provided." );
		utility::vector1< string > const codes_files( 1, "sentence_case.codes" );
		option[ OptionKeys::in::alternate_3_letter_codes ]( codes_files );

		residue_empty = NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "" );
		residue_sentence_case = NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "Ala" );
		residue_all_caps = NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "ALA" );
		residue_special = NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "Hcy" );
		residue_bogus = NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "WOO!" );

		TS_ASSERT_EQUALS( residue_empty.first, "" );
		TS_ASSERT_EQUALS( residue_empty.second, "" );
		TS_ASSERT_EQUALS( residue_sentence_case.first, "ALA" );
		TS_ASSERT_EQUALS( residue_sentence_case.second, "" );
		TS_ASSERT_EQUALS( residue_all_caps.first, "ALA" );
		TS_ASSERT_EQUALS( residue_all_caps.second, "" );
		TS_ASSERT_EQUALS( residue_special.first, "HCY" );
		TS_ASSERT_EQUALS( residue_special.second, "HOMOCYSTEINE" );
		TS_ASSERT_EQUALS( residue_bogus.first, "WOO!" );
		TS_ASSERT_EQUALS( residue_bogus.second, "" );
	}
};  // class NomenclatureManagerTests
