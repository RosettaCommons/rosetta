// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/utility.cxxtest.hh
/// @brief  utility.cxxtest: placeholder for main() function when creating unit test executable
/// @author Sergey Lyskov


#ifndef INCLUDED_utility_utility_cxxtest_HH
#define INCLUDED_utility_utility_cxxtest_HH

// When compiling tests for CxxTest using the --part argument to the Python script, no main()
// method is added to the generated .cpp files. But there needs to be one file that has a main.
// Let's assume inside each "project" there will be a project.cxxtest.hh which has no tests or anything
// in it. This file will get processed with the --root switch which adds a main() method to the
// generated .cpp file.  This is that file for the utility project.

#include <cxxtest/TestSuite.h>
#include <utility/exit.hh>


//class UtilityTests : public CxxTest::TestSuite {
//public:

// void test_utility_exit() {
// 	utility_exit_with_message("Testing utility exit... this test should fail...");
// }

//	void test_assert_exit() {
//		// Test to see if a failed assert is registered as a fail on the test server
//		assert( 0 );
//	}

//};


#endif //INCLUDED_utility_cxxtest_HH
