// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/file/file_sys_util.cxxtest.hh
/// @brief  file_sys_util unit test suite
/// @author Ian Davis

// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/file/file_sys_util.hh>
#include <utility/file/PathName.hh>

// C++ headers
#include <iostream>

class FileSysUtilTests : public CxxTest::TestSuite {

	public:

	void test_create_directory_recursive() {
		using namespace utility::file;
		PathName p1("foo/bar/baz");
		TS_ASSERT( create_directory_recursive(p1.name()) );
		TS_ASSERT( file_exists(p1.name()) );
	}

};

