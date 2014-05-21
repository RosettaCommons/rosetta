// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/numbers.cxxtest.hh
/// @brief  test suite for simple utilities on numbers
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// Package headers
#include <utility/numbers.hh>

class UtilityNumbersTests : public CxxTest::TestSuite {

public:
	void setUp() {
	}

	void test_undefined() {
		TS_ASSERT( utility::is_undefined( utility::get_undefined_size() ) );
		TS_ASSERT( utility::is_undefined( utility::get_undefined_real() ) );
		TS_ASSERT( ! utility::is_undefined( platform::Size(0) ) );
		TS_ASSERT( ! utility::is_undefined( 0.0 ) );
	}

}; // NumericUtilTests
