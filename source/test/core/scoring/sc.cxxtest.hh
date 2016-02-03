// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/sc.cxxtest.hh
/// @brief  test suite for core::scoring::sc
/// @author Luki Goldschmidt

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/pose/Pose.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>

// Package Headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;

class ScTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}


	// Shared finalization goes here.
	void tearDown() {}


	// --------------- Test Cases --------------- //
	void test_sc() {

		core::scoring::sc::ShapeComplementarityCalculator scc;

		struct {
			const char *fn;
			core::Size atoms_0;
			core::Size atoms_1;
			float sc;
			float sc_quick;
			float area;
			float distance;
		} test_structures[] = {
		{ "core/scoring/sc_NNQQNY.pdb", 271, 271, 0.714775, 0.720340, 950.14, 0.5008 }
			};

		for ( unsigned int i = 0;
				i < sizeof(test_structures)/sizeof(*test_structures);
				++i ) {
			core::pose::Pose pose;
			int r;

			core::import_pose::pose_from_file( pose, test_structures[i].fn , core::import_pose::PDB_file);

			// Test standard mode
			{
				scc.Reset();
				// scc.settings.density = 5; // Quick mode

				r = scc.Calc(pose, 1);
				TS_ASSERT_EQUALS(r, 1);

				TS_ASSERT_EQUALS( scc.GetResults().surface[0].nAtoms, test_structures[i].atoms_0 );
				TS_ASSERT_EQUALS( scc.GetResults().surface[1].nAtoms, test_structures[i].atoms_1 );

				TS_ASSERT_DELTA( scc.GetResults().sc, test_structures[i].sc, 0.001 );
				TS_ASSERT_DELTA( scc.GetResults().area, test_structures[i].area, 1 );
				TS_ASSERT_DELTA( scc.GetResults().distance, test_structures[i].distance, 0.001 );
			}

			// Test simple calculation via static method
			{
				scc.Reset();
				// scc.settings.density = 5; // Quick mode

				core::Real sc = core::scoring::sc::ShapeComplementarityCalculator::CalcSc(pose, 1, 1);
				TS_ASSERT_DELTA( sc, test_structures[i].sc_quick, 0.001 );
			}
		}
	}
};


