// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cyclic_peptide/AchiralResidueTypeTests.cxxtest.hh
/// @brief  Unit tests for achiral residue types other than glycine (e.g. 2-aminoisobutyric acid [AIB]).
/// @details This version tests using the default score function (currently talaris2014).
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <test/protocols/cyclic_peptide/AchiralResidueTypeTestHeaders.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("AchiralResidueTypeTests");


class AchiralResidueTypeTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-symmetric_gly_tables true" );
	}

	void tearDown(){
	}

	void test_mirror_symmetry_AIB(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_mirror_symmetry( TR, "AIB", false );
	}

	void test_mirror_symmetry_AIB_before_proline(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_mirror_symmetry( TR, "AIB", true );
	}

	void test_mirror_symmetry_GLY(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_mirror_symmetry( TR, "GLY", false );
	}

	void test_mirror_symmetry_GLY_before_proline(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_mirror_symmetry( TR, "GLY", true );
	}

	void test_symmetric_rama_prepro_scoring_AIB(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_symmetric_rama_prepro_scoring( TR, "AIB", false );
	}

	void test_symmetric_rama_prepro_scoring_AIB_before_proline(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_symmetric_rama_prepro_scoring( TR, "AIB", true );
	}

	void test_symmetric_rama_prepro_scoring_GLY(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_symmetric_rama_prepro_scoring( TR, "GLY", false );
	}

	void test_symmetric_rama_prepro_scoring_GLY_before_proline(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_symmetric_rama_prepro_scoring( TR, "GLY", true );
	}

};



