// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  test/core/io/pdb/NomenclatureManagerTests.cxxtest.hh
/// @brief   Test suite for the NomenclatureManager singleton
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/io/NomenclatureManager.hh>

// Utility header
#include <utility/vector1.hh>

// Basic headers
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

//read pose stuff
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/pdb/pdb_reader.hh>
#include <core/io/StructFileReaderOptions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// C++ header
#include <string>

static THREAD_LOCAL basic::Tracer TR("core.io.NomenclatureManagerTests.cxxtest");

// Without the command line (wocl)
class NomenclatureManagerTestsWOCL : public CxxTest::TestSuite {
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
		using namespace core::io;

		TR <<  "Testing rosetta_names_from_pdb_code() static method with alternative 3-letter codes not provided."  << std::endl;

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
	}
};

// With the command line (wcl)
class NomenclatureManagerTestsWCL : public CxxTest::TestSuite {
public: // Standard methods ///////////////////////////////////////////////////
	void setUp()
	{
		core_init();
	}

	void test_nomenclature_manager_test_w_command_line_alt_codes() {

		using namespace std;
		using namespace basic::options;
		using namespace core::io;


		utility::vector1< string >  codes_files( 1, "sentence_case.codes" );
		option[ OptionKeys::in::alternate_3_letter_codes ]( codes_files );

		typedef pair< string, string > string_pair;

		string_pair residue_empty = NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "" );
		string_pair residue_sentence_case = NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "Ala" );
		string_pair residue_all_caps = NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "ALA" );
		string_pair residue_special = NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "Hcy" );
		string_pair residue_bogus = NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "WOO!" );

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

	void test_nomenclature_manager_three_name_conflicts(){
		using namespace std;
		using namespace basic::options;
		using namespace core::io;

		TR <<  "Testing nomenclature_manager_three_name_conflicts() static method with alternative 3-letter codes provided."  << std::endl;
		//        utility::vector1< string > const codes_files( 1, "gylcam.codes" );
		//        option[ OptionKeys::in::alternate_3_letter_codes ]( codes_files );
		//
		//        typedef pair< string, string > string_pair;
		//
		//        string_pair onega_to_glc =  NomenclatureManager::get_instance()->rosetta_names_from_pdb_code( "1GA" );
		//
		//        //make sure that this three letter name is the same for the glucose
		//        TS_ASSERT_EQUALS( onega_to_glc.second, "Glc");
		//
		//
		//
		//        core::pose::Pose pose;
		//        core::io::StructFileReaderOptions options;
		//        core::io::StructFileRep sfr = core::io::pdb::create_sfr_from_pdb_file_contents( "core/io/1l2y_w_glucose.pdb", options );
		//
		//        chemical::ResidueTypeSetCOP residue_set
		//        ( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );
		//        PoseFromSFRBuilder pb( residue_set, options );
		//        pb.build_pose( sfr, pose );
		//
		//        TS_ASSERT( pose.total_residue() == 20 );
		//        if ( pose.total_residue() != 20 ) return;
		//
		//        TS_ASSERT( pose.residue_type(  1 ).has_variant_type( chemical::LOWER_TERMINUS_VARIANT ));
		//        for ( Size ii = 2; ii <= 19; ++ii ) {
		//            TS_ASSERT( ! pose.residue_type( ii ).has_variant_type( chemical::LOWER_TERMINUS_VARIANT ));
		//            TS_ASSERT( ! pose.residue_type( ii ).has_variant_type( chemical::UPPER_TERMINUS_VARIANT ));
		//        }
		//        TS_ASSERT( pose.residue_type( 20 ).has_variant_type( chemical::UPPER_TERMINUS_VARIANT ));


	}

};  // class NomenclatureManagerTests
