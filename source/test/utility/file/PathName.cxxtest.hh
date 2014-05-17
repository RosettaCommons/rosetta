// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/file/PathName.cxxtest.hh
/// @brief  PathName unit test suite
/// @author Ian Davis

// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/file/PathName.hh>

// C++ headers
#include <iostream>

class PathNameTests : public CxxTest::TestSuite {

	public:

	/// @brief Test function for determining parent directory of a given path
	void test_parent() {
		using namespace utility::file;
		// Absolute paths
		PathName p1("/foo/bar/baz");
		//std::cout << p1.name() << std::endl;
		//std::cout << p1.parent().name() << std::endl;
		//std::cout << p1.parent().parent().name() << std::endl;
		//std::cout << p1.parent().parent().parent().name() << std::endl;
		//std::cout << p1.parent().parent().parent().parent().name() << std::endl;
		TS_ASSERT_EQUALS( PathName("/foo/bar"), p1.parent() );
		TS_ASSERT_EQUALS( PathName("/foo"),     p1.parent().parent() );
		TS_ASSERT_EQUALS( PathName("/"),        p1.parent().parent().parent() );
		TS_ASSERT_EQUALS( PathName("/"),        p1.parent().parent().parent().parent() );
		// Relative paths
		PathName p2("foo/bar/baz");
		//std::cout << p2.name() << std::endl;
		//std::cout << p2.parent().name() << std::endl;
		//std::cout << p2.parent().parent().name() << std::endl;
		//std::cout << p2.parent().parent().parent().name() << std::endl;
		//std::cout << p2.parent().parent().parent().parent().name() << std::endl;
		TS_ASSERT_EQUALS( PathName("foo/bar"), p2.parent() );
		TS_ASSERT_EQUALS( PathName("foo"),     p2.parent().parent() );
		TS_ASSERT_EQUALS( PathName(""),        p2.parent().parent().parent() );
		TS_ASSERT_EQUALS( PathName(""),        p2.parent().parent().parent().parent() );
	}

	/// @brief Test function for path assignment
	void test_path() {
		using namespace utility::file;
		// Case 1 paths
		PathName pa1("/foo/bar");
		TS_ASSERT_EQUALS( PathName("\\foo\\bar").path(), pa1.path() );
		// Case 2 paths
		PathName pa2("\\foo\\bar");
		TS_ASSERT_EQUALS( PathName("/foo/bar").path(), pa2.path() );
	}
	/// @brief Test function for volume assigment
	void test_vol() {
		using namespace utility::file;
		// Case 1 with volume stripping
		PathName pa3("C:/foo/bar");
		TS_ASSERT_EQUALS( PathName("\\foo\\bar").vol(), pa3.vol() );
		// Case 2 with volume stripping
		PathName pa4("C:\\foo\\bar");
		TS_ASSERT_EQUALS( PathName("\\foo\\bar").vol(), pa4.vol() );
	}
};

