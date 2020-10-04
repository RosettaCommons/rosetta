// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    basic/execute.cxxtest.hh
///
/// @brief   unit test for basic::execute function
///
/// @author  Sergey Lyskov

// Test headers
#include <cxxtest/TestSuite.h>

#include <basic/execute.hh>

class ExecuteTest : public CxxTest::TestSuite {
public:
	void setUp() {}

	void tearDown() {}

	void test_shell() {
		using namespace basic;

		basic::ExecutionResult r;
		r = execute("Testing ls...", "ls", {}, false, true);
		TS_ASSERT_EQUALS(r.result, 0);

		r = execute("Testing ls...", "ls && cd .", {}, false, true);
		TS_ASSERT_DIFFERS(r.result, 0);
		TS_ASSERT_STRING_CONTAINS(r.output, "ERROR: execvp failed");
	}
};
