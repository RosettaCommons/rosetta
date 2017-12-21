// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  utility/json/json_utilitiesTests.cxxtest.hh
/// @brief  test nlohmann::json utilities.  test var references are to "jenny" because if I am asked to think of a random number for a unit test I invariably think of 867-5309
/// @author Steven Lewis (smlewi@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers

// Core Headers

// Utility, etc Headers
#include <utility/json/json_utilities.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("test.utility.json.json_utilitiesTests");

using namespace utility::json;

class json_utilitiesTests : public CxxTest::TestSuite {

public:
	std::string const thingname = "thingname";
	std::string const missingname = "missingname";
	std::string const missing_message ="JSON element missingname missing";
	std::string const wrongtype_message ="JSON element thingname not of ";

	void setUp(){
		core_init();
	}

	void tearDown(){
	}

	void test_bool(){

		json const test_json({{thingname, true}});

		//TR << test_json << std::endl;

		bool bool_check(false);
		extract_boolean_from_json(test_json, thingname, bool_check); //check that extraction succeeds

		TS_ASSERT_EQUALS(bool_check, true);

		try { //check that missing extraction fails in the right way
			set_throw_on_next_assertion_failure();
			extract_boolean_from_json(test_json, missingname, bool_check);
			TS_ASSERT(false); //should have thrown
		} catch ( utility::excn::Exception const & e ) {
			std::string const & message(e.msg());
			//std::cout << "message: " << message << std::endl;
			TS_ASSERT_DIFFERS(message.find(missing_message), std::string::npos); //check that the exception message contains the expected string
		}

		json const wrongtype_json({{thingname, "jenny"}});
		try { //check that wrong-type extraction fails in the right way.
			set_throw_on_next_assertion_failure();
			extract_boolean_from_json(wrongtype_json, thingname, bool_check);
			TS_ASSERT(false); //should have thrown
		} catch ( utility::excn::Exception const & e ) {
			std::string const & message(e.msg());
			//std::cout << "message: " << message << std::endl;
			TS_ASSERT_DIFFERS(message.find(wrongtype_message), std::string::npos); //check that the exception message contains the expected string
		}

		//sanity check: does this library eliminate duplicate entries? (yes)
		//json const test2_json({{thingname, true}, {thingname, false}});
		//TR << test2_json << std::endl;
		//extract_boolean_from_json(test2_json, thingname, bool_check);

	};

	void test_number(){

		platform::Real const digits(867.5309);
		json const test_json({{thingname, digits}});

		TR << test_json << std::endl;

		platform::Real number_check(0.0);
		extract_number_from_json(test_json, thingname, number_check); //check that extraction succeeds

		TS_ASSERT_EQUALS(number_check, digits); //should be an ok float comparison; delta it (or use a precisely representable, non-song-referencing number) as needed

		try { //check that missing extraction fails in the right way
			set_throw_on_next_assertion_failure();
			extract_number_from_json(test_json, missingname, number_check);
			TS_ASSERT(false); //should have thrown
		} catch ( utility::excn::Exception const & e ) {
			std::string const & message(e.msg());
			//std::cout << "message: " << message << std::endl;
			TS_ASSERT_DIFFERS(message.find(missing_message), std::string::npos); //check that the exception message contains the expected string
		}

		json const wrongtype_json({{thingname, "jenny"}});
		try { //check that wrong-type extraction fails in the right way.
			set_throw_on_next_assertion_failure();
			extract_number_from_json(wrongtype_json, thingname, number_check);
			TS_ASSERT(false); //should have thrown
		} catch ( utility::excn::Exception const & e ) {
			std::string const & message(e.msg());
			//std::cout << "message: " << message << std::endl;
			TS_ASSERT_DIFFERS(message.find(wrongtype_message), std::string::npos); //check that the exception message contains the expected string
		}

		//also check extraction of core::Size types
		platform::Size const int_digits(867);
		extract_number_from_json(test_json, thingname, number_check);  //still have to extract as Real
		platform::Size const int_number_check(number_check); // then cast to Size
		TS_ASSERT_EQUALS(int_digits, int_number_check);

	};

	void test_nonempty_string(){

		std::string const digits("867.5309");
		json const test_json({{thingname, digits}});

		TR << test_json << std::endl;

		std::string nonempty_string_check("jenny");
		extract_nonempty_string_from_json(test_json, thingname, nonempty_string_check); //check that extraction succeeds

		TS_ASSERT_EQUALS(nonempty_string_check, digits);

		try { //check that missing extraction fails in the right way
			set_throw_on_next_assertion_failure();
			extract_nonempty_string_from_json(test_json, missingname, nonempty_string_check);
			TS_ASSERT(false); //should have thrown
		} catch ( utility::excn::Exception const & e ) {
			std::string const & message(e.msg());
			//std::cout << "message: " << message << std::endl;
			TS_ASSERT_DIFFERS(message.find(missing_message), std::string::npos); //check that the exception message contains the expected string
		}

		json const wrongtype_json({{thingname, false}});
		try { //check that wrong-type extraction fails in the right way.
			set_throw_on_next_assertion_failure();
			extract_nonempty_string_from_json(wrongtype_json, thingname, nonempty_string_check);
			TS_ASSERT(false); //should have thrown
		} catch ( utility::excn::Exception const & e ) {
			std::string const & message(e.msg());
			//std::cout << "message: " << message << std::endl;
			TS_ASSERT_DIFFERS(message.find(wrongtype_message), std::string::npos); //check that the exception message contains the expected string
		}

	};

	void test_nonempty_array(){

		json const author_array({"by", "tommy", "tutone"});
		json const test_json({{thingname, author_array}});

		TR << test_json << std::endl;

		json nonempty_array_check;
		extract_nonempty_array_from_json(test_json, thingname, nonempty_array_check); //check that extraction succeeds

		TS_ASSERT_EQUALS(nonempty_array_check, author_array);

		try { //check that missing extraction fails in the right way
			set_throw_on_next_assertion_failure();
			extract_nonempty_array_from_json(test_json, missingname, nonempty_array_check);
			TS_ASSERT(false); //should have thrown
		} catch ( utility::excn::Exception const & e ) {
			std::string const & message(e.msg());
			//std::cout << "message: " << message << std::endl;
			TS_ASSERT_DIFFERS(message.find(missing_message), std::string::npos); //check that the exception message contains the expected string
		}

		json const wrongtype_json({{thingname, false}});
		try { //check that wrong-type extraction fails in the right way.
			set_throw_on_next_assertion_failure();
			extract_nonempty_array_from_json(wrongtype_json, thingname, nonempty_array_check);
			TS_ASSERT(false); //should have thrown
		} catch ( utility::excn::Exception const & e ) {
			std::string const & message(e.msg());
			//std::cout << "message: " << message << std::endl;
			TS_ASSERT_DIFFERS(message.find(wrongtype_message), std::string::npos); //check that the exception message contains the expected string
		}

	};

	void test_nonempty_object(){

		json const song_object({{"by", "tommy tutone"}, {"number", 867.5309}});
		json const test_json({{thingname, song_object}});

		TR << test_json << std::endl;

		json nonempty_object_check;
		extract_nonempty_object_from_json(test_json, thingname, nonempty_object_check); //check that extraction succeeds

		TS_ASSERT_EQUALS(nonempty_object_check, song_object);

		try { //check that missing extraction fails in the right way
			set_throw_on_next_assertion_failure();
			extract_nonempty_object_from_json(test_json, missingname, nonempty_object_check);
			TS_ASSERT(false); //should have thrown
		} catch ( utility::excn::Exception const & e ) {
			std::string const & message(e.msg());
			//std::cout << "message: " << message << std::endl;
			TS_ASSERT_DIFFERS(message.find(missing_message), std::string::npos); //check that the exception message contains the expected string
		}

		json const wrongtype_json({{thingname, false}});
		try { //check that wrong-type extraction fails in the right way.
			set_throw_on_next_assertion_failure();
			extract_nonempty_object_from_json(wrongtype_json, thingname, nonempty_object_check);
			TS_ASSERT(false); //should have thrown
		} catch ( utility::excn::Exception const & e ) {
			std::string const & message(e.msg());
			//std::cout << "message: " << message << std::endl;
			TS_ASSERT_DIFFERS(message.find(wrongtype_message), std::string::npos); //check that the exception message contains the expected string
		}

	};


};
