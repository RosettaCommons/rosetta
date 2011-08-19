// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/string_util.cxxtest.hh
/// @brief  string_util.cxxtest: test suite for utility::string_util
/// @author James Thompson
/// @author Christopher Miles (cmiles@uw.edu)

// Testing headers
#include <cxxtest/TestSuite.h>

// Project headers
#include <utility/string_util.hh>

// C/C++ headers
#include <iostream>
#include <string>
#include <sstream>

using std::endl;
using std::string;
using std::stringstream;

class StringUtilTests : public CxxTest::TestSuite {
 public:

	// duplicated implementation in FileName... now it's just an alias for utility::basename
	void test_file_basename() {
		TS_ASSERT_EQUALS(utility::file_basename("core/scoring/ScoreFunction.cc"),
						 "ScoreFunction.cc");
	}

	void test_filename() {
		TS_ASSERT_EQUALS(utility::filename("/foo/bar/baz"), "baz");
	}

	void test_pathname() {
		TS_ASSERT_EQUALS(utility::pathname("/foo/bar/baz"), "/foo/bar/");
	}

	void test_same_ignoring_spaces() {
		TS_ASSERT(utility::same_ignoring_spaces("CA", "CA"));
		TS_ASSERT(utility::same_ignoring_spaces(" CA", "CA"));
		TS_ASSERT(utility::same_ignoring_spaces("CA ", "CA"));
		TS_ASSERT(utility::same_ignoring_spaces("   CA     ", "  CA  "));
	}

	void test_string2uint() {
		string s = "42";
		unsigned int i;
		utility::string2uint(s, &i);
		TS_ASSERT_EQUALS(i, 42);
	}

	void test_readfromfileordie() {
		string contents;
		utility::ReadFromFileOrDie("utility/io/no_final_newline.txt", &contents);

		stringstream result;
		result << "foo" << endl
		       << "bar" << endl
		       << "baz has no newline!" << endl;

		TS_ASSERT_EQUALS(contents, result.str());
	}
};
