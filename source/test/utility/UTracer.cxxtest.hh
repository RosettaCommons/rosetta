// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/UTracer.cxxtest.hh
/// @brief  Test for the UTracer
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Package headers
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// This really isn't a test of Rosetta, but rather a test to make sure that the UTracer unit testing hook works
// how expected. (It's only in the utility namespace because there isn't specifically a unit testing subtest.

class UTracerTests : public CxxTest::TestSuite {

	public:

		void test_comparisions() {
			test::UTracer UT("utility/UTracer_test.u");

			UT << "These lines will be compared exactly" << std::endl;
			UT << "Spaces 21 9 1.7 9e-16 0.2 3.14159 -0.726335 -2.43169 -0.834851 0 -1.75902 3.90664e-16" << std::endl;
			UT << "Tabs 21\t9\t1.7\t9e-16\t0.2\t3.14159\t-0.726335\t-2.43169\t-0.834851\t0\t-1.75902\t3.90664e-16" << std::endl;
			UT.abs_tolerance(1e-15);
			UT << std::endl;
			UT << "These lines will be compared approximately" << std::endl;
			UT << "NOTE: The discrepancy between the in-code and file representations is *deliberate* - please don't overwrite willy-nilly." << std::endl;
			// In code this is the same as the lines above. In the comparison file, we've altered the values at just past the cutoff limit.
			UT << "Spaces 21 9 1.7 9e-16 0.2 3.14159 -0.726335 -2.43169 -0.834851 0 -1.75902 3.90664e-16" << std::endl;
			UT << "Tabs 21\t9\t1.7\t9e-16\t0.2\t3.14159\t-0.726335\t-2.43169\t-0.834851\t0\t-1.75902\t3.90664e-16" << std::endl;
		}

};


