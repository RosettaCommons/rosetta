// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/file/FileName.cxxtest.hh
/// @brief  FileName unit test suite
/// @author Amanda Duran Alyssa Lokits

// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/file/FileName.hh>

// C++ headers
#include <iostream>

class FileNameTests : public CxxTest::TestSuite {

	public:

	/// @brief Test function for path assignment
	void test_path() {
		using namespace utility::file;
		// Case 1 paths
		FileName f1("/foo/bar");
		TS_ASSERT_EQUALS( FileName("\\foo\\bar").path(), f1.path() );
		// Case 2 paths
		FileName f2("\\foo\\bar");
		TS_ASSERT_EQUALS( FileName("/foo/bar").path(), f2.path() );
	}
	/// @brief Test function for volume assigment
	void test_vol() {
		using namespace utility::file;
		// Case 1 with volume stripping
		FileName v1("C:/foo/bar");
		TS_ASSERT_EQUALS( FileName("\\foo\\bar").vol(), v1.vol() );
		// Case 2 with volume stripping
		FileName v2("C:\\foo\\bar");
		TS_ASSERT_EQUALS( FileName("\\foo\\bar").vol(), v1.vol() );
	}

};

