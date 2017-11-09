// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /rosetta/main/source/test/protocols/legacy_sewing/SewingHasherTests.cxxtest.hh
/// @brief test suite for sewing_hasher
/// @author Doonam Kim


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit headers
#include <protocols/legacy_sewing/hashing/Hasher.hh>
#include <test/core/init_util.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/sql_database/DatabaseSessionManager.hh>

/// Project headers

//Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

//Devel headers
#include <devel/init.hh>

// --------------- Test Class --------------- //
using namespace protocols::legacy_sewing;

static basic::Tracer tr("protocols.legacy_sewing.SewingHasherTests.cxxtest");

class SewingHasherTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		core_init();

	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_key_generation() {
		//devel::sewing::SewHash hasher;
		//TS_ASSERT_EQUALS(false, hasher.hash_map().key_eq()(key2, key3));
	}

	void test_transform_coords() {
		//devel::sewing::SewHash hasher;
	}

	void test_scoring() {
	}

	void test_serialization() {
		//TS_ASSERT_EQUALS(hasher.hash_map().size(), hasher2.hash_map().size());
	}

	void test_5_ss_model_generation_from_db() {
		// tr << "comments should work" << std::endl;

		using utility::sql_database::DatabaseSessionManager;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		//Create comments stream for model file and add the date
		std::stringstream comments;
		time_t t = time(0);   // get time now
		struct tm * now = localtime( & t );
		comments << "#Model file created on " << (now->tm_year + 1900) << '-'
			<< (now->tm_mon + 1) << '-'
			<<  now->tm_mday
			<< std::endl;

		//Generate models from a features database. Each segment is a single piece of secondary structure
		option[inout::dbms::database_name].value("protocols/legacy_sewing/inputs/1ten_m8.db3");

		tr << "#Model will be generated from sqlite database " << option[inout::dbms::database_name].value() << std::endl;

		//Generate models from a features database. Each segment is a single piece of secondary structure
		comments << "#Models generated from sqlite database " << option[inout::dbms::database_name].value() << std::endl;
		std::string hash_between = "hash_tag_terminal_HEs";
		std::map<int, Model> models = get_5_ss_models_from_db(hash_between); //may work for now

		std::string model_five_ss_filename = "unit_test_model_w_five_ss";
		write_model_file(comments.str(), models, model_five_ss_filename);
		tr << "New model file with 5 ss " << model_five_ss_filename << " successfully written." << std::endl;

		std::ifstream file("unit_test_model_w_five_ss");
		std::string str;
		//    std::string file_contents;
		int line_i=0;
		while ( std::getline(file, str) ) {
			line_i++;
			//tr << "line_i: " << line_i << std::endl;
			//tr << "str: " << str << std::endl;
			if ( line_i == 895 ) {
				tr << "line_i: " << line_i << std::endl;
				tr << "str: " << str << std::endl;
				TS_ASSERT_EQUALS(str, "ATOM 4 50.997 48.694 37.636");
			} else if ( line_i == 896 ) {
				tr << "line_i: " << line_i << std::endl;
				tr << "str: " << str << std::endl;
				TS_ASSERT_EQUALS(str, "RESIDUE 63 0 ALA ");
			} else if ( line_i == 897 ) {
				tr << "line_i: " << line_i << std::endl;
				tr << "str: " << str << std::endl;
				TS_ASSERT_EQUALS(str, "ATOM 1 49.093 49.4 38.577");
			}
		}
	}//test_5_ss_model_generation_from_db

private:

};
