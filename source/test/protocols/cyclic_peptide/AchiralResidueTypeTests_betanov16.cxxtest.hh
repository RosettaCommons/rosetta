// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cyclic_peptide/AchiralResidueTypeTests_betanov16.cxxtest.hh
/// @brief  Unit tests for achiral residue types other than glycine (e.g. 2-aminoisobutyric acid [AIB]).
/// @details This version tests using the beta_nov16 score function.
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

static basic::Tracer TR("AchiralResidueTypeTests_betanov16");


class AchiralResidueTypeTests_betanov16 : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-beta_nov16 -symmetric_gly_tables true" );
	}

	void tearDown(){
	}

	void test_import_AIB(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_import( TR, "protocols/cyclic_peptide/AIB_pose.pdb", 10, "AIB", "AIB" );
	}

	void test_mirror_symmetry_AIB(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_mirror_symmetry( TR, "AIB", false, false );
	}

	void test_mirror_symmetry_AIB_before_proline(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_mirror_symmetry( TR, "AIB", true, false );
	}

	void test_mirror_symmetry_GLY(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_mirror_symmetry( TR, "GLY", false, false );
	}

	void test_mirror_symmetry_GLY_before_proline(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_mirror_symmetry( TR, "GLY", true, false );
	}

	void test_mirror_symmetry_sarcosine(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_mirror_symmetry( TR, "GLY", false, true );
	}

	void test_mirror_symmetry_sarcosine_before_proline(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_mirror_symmetry( TR, "GLY", true, true );
	}

	void test_symmetric_rama_prepro_scoring_AIB(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_symmetric_rama_prepro_scoring( TR, "AIB", false, false, 0.000001 );
	}

	void test_symmetric_rama_prepro_scoring_AIB_before_proline(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_symmetric_rama_prepro_scoring( TR, "AIB", true, false, 0.000001 );
	}

	void test_symmetric_rama_prepro_scoring_GLY(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_symmetric_rama_prepro_scoring( TR, "GLY", false, false, 0.000001 );
	}

	void test_symmetric_rama_prepro_scoring_GLY_before_proline(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_symmetric_rama_prepro_scoring( TR, "GLY", true, false, 0.000001);
	}

	void test_symmetric_rama_prepro_scoring_sarcosine(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_symmetric_rama_prepro_scoring( TR, "GLY", false, true, 0.000001 );
	}

	void test_symmetric_rama_prepro_scoring_sarcosine_before_proline(){
		test::protocols::cyclic_peptide::AchiralResidueTypeTestHelper helper;
		helper.test_symmetric_rama_prepro_scoring( TR, "GLY", true, true, 0.000001 );
	}

};



