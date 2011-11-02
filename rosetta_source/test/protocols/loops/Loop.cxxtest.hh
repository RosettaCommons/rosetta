// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loops/Loop.cxxtest.hh
/// @brief test suite for protocols/loops/Loop
/// @author Christopher Miles (cmiles@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

// Project headers
#include <protocols/loops/Loop.hh>

//Auto Headers
#include <utility/vector1.hh>


namespace {

class LoopTest : public CxxTest::TestSuite {
 public:
	void test_length() {
		// default constructor
		protocols::loops::Loop l1;
		TS_ASSERT_EQUALS(l1.length(), 1);

		// input constructor
		protocols::loops::Loop l2(3, 8);
		TS_ASSERT_EQUALS(l2.length(), 6);
	}
};
}  // anonymous namespace
