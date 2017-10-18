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
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/core/scoring/RamaPrePro_util.h>

// Project Headers
#include <core/scoring/RamaPrePro.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ResidueTypeSet.hh>

// Protocol Headers
#include <protocols/cyclic_peptide/FlipChiralityMover.hh>

// Basic Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("RamaPreProTests_symmgly");


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



