// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file utility/options/OptionCollection.cxxtest.hh
/// @brief test suite for options system
/// @author Matthew O'Meara mattjomeara@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>

class OptionsSystemTests : public CxxTest::TestSuite {

public:

	void setUp() {
		basic::options::initialize();
	}

	void test_find_key_cl_good(){
		using namespace basic::options;
		try{
			std::string key(option.find_key_cl("in:file:s", "", true));
		} catch (...){
			TS_ASSERT(false);
		}

		try{
			std::string key(option.find_key_cl("s", "in:file", false));
		} catch (...){
			TS_ASSERT(false);
		}

		try{
			std::string key(option.find_key_cl("s", "", false));
		} catch (...){
			TS_ASSERT(false);
		}

		try{
			std::string key(option.find_key_cl("s", "", true));
		} catch (...){
			TS_ASSERT(false);
		}


		try{
			std::string key(option.find_key_cl("in:file:s", "", false));
		} catch (...){
			TS_ASSERT(false);
		}


	}

	void test_find_key_cl_error() {
		using namespace basic::options;
		try{
			std::string key(option.find_key_cl(":::s", "in:file", true));
		} catch (...){
			TS_ASSERT(true);
		}

		try{
			// unique best suffix match
			std::string key(option.find_key_cl("asdfasdfasf:s", "in:file", false));
		} catch (...){
			TS_ASSERT(true);
		}

		try{
			// ignore context when 'top' == true
			std::string key(option.find_key_cl("s", "in:file", true));
		} catch (...){
			TS_ASSERT(true);
		}

		try{
			std::string key(option.find_key_cl("in::file::s", "", false));
		} catch (...){
			TS_ASSERT(true);
		}

		try{
			std::string key(option.find_key_cl("-in::file::s-", "", false));
		} catch (...){
			TS_ASSERT(true);
		}

		try{
			std::string key(option.find_key_cl("in:file", "", false));
		} catch (...){
			TS_ASSERT(true);
		}

		try{
			std::string key(option.find_key_cl("-in:file:s", "", false));
		} catch (...){
			TS_ASSERT(true);
		}

	}
};
