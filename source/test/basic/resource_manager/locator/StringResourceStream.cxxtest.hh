// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/StringResourceStream.cxxtest.hh
/// @brief  Test the StringResourceStream class
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test Headers
#include <test/core/init_util.hh>
#include <cxxtest/TestSuite.h>
#include <basic/Tracer.hh>

#include <basic/resource_manager/locator/StringResourceStream.hh>

// Utility headers

// C++ headers
#include <sstream>

static basic::Tracer tr( "basic.resource_manager.locator.StringResourceStream.cxxtest");

class StringResourceStreamTests : public CxxTest::TestSuite {

public:

	void test_null() {
	}

	//  void test_default_construction() {
	//
	//  using namespace basic::resource_manager;
	//  using namespace basic::resource_manager::locator;
	//
	//    StringResourceStream srs;
	//    srs.fill("abc");
	//
	//    std::stringstream out;
	//    out << srs.stream();
	//
	//    std::string out_str(out.str());
	//    TS_ASSERT_EQUALS("abc", out_str);
	//  }

	// TODO: MAKE THIS TEST WORK
	//
	//  void test_stream_construction() {
	//
	//  using namespace basic::resource_manager;
	//  using namespace basic::resource_manager::locator;
	//
	//    std::stringstream in;
	//    in << "abc";
	//
	//    StringResourceStream srs(in);
	//
	//    std::stringstream out;
	//    out << srs.stream();
	//
	//    TS_ASSERT_EQUALS("abc", out.str());
	//  }

};
