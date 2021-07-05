// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/RamaPrePro.cxxtest.hh
/// @brief  Unit tests for the RamaPrePro energy.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/scoring/RamaPrePro_util.hh>

// Project Headers

// Core Headers

// Protocol Headers

// Basic Headers
#include <basic/Tracer.hh>

#include <core/init_util.hh> // AUTO IWYU For core_init_with_additional_options

static basic::Tracer TR("RamaPreProTests_symmgly");


class RamaPreProTests_symmgly : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options("-symmetric_gly_tables true");
	}

	void tearDown(){
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution of gly.
	void test_random_phipsi_gly() {
		do_ramaprepro_test("AGAA", false, false);
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution for a pre-proline gly.
	void test_random_phipsi_gly_prepro() {
		do_ramaprepro_test("AGPA", false, false);
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution of gly, with flipping.
	void test_random_phipsi_gly_invert() {
		do_ramaprepro_test("AGAA", false, true);
	}

	/// @brief Test the drawing of random mainchain torsion values from the
	/// Ramachandran probability distribution for a pre-proline gly, with flipping.
	void test_random_phipsi_gly_prepro_invert() {
		do_ramaprepro_test("AGPA", false, true);
	}

};



