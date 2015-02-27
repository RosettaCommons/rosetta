// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/P_AA_Energy.cxxtest.hh
/// @brief  test suite for core::scoring::P_AA_Energy.cc
/// @author Ron Jacak

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/methods/P_AA_Energy.hh>

#include <platform/types.hh>

// Package Headers
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <basic/options/option.hh>

//Auto Headers
#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>



// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class P_AA_EnergyTests : public CxxTest::TestSuite {

	public:

	PoseOP the_pose;
	P_AA_EnergyOP paa_energy;


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {

		using namespace std;

		using namespace core;
		core_init();

		the_pose = create_test_in_pdb_poseop();
		//core::import_pose::pose_from_pdb( *the_pose, "core/scoring/methods/test_in.pdb" );

		paa_energy = P_AA_EnergyOP( new P_AA_Energy );

	}

	// Shared finalization goes here.
	void tearDown() {
		the_pose.reset();
		paa_energy.reset();
	}


	// --------------- Test Cases --------------- //
	void test_eval_energy() {

		float const TOLERATED_ERROR = 0.0001;

		Real correct_answers[] = {
			2.8283, 2.4895, 2.8936, 2.8335, 2.8936, 3.8162, 2.7905, 2.8936, 2.4973, 2.8283, 4.2184, 2.8936, 2.7931,
			2.8283, 3.0597, 2.4973, 2.7931, 2.7905, 3.0627, 2.4973, 2.7905, 2.4973, 2.7931, 2.8163, 2.6607, 2.7905,
			2.7931, 3.0664, 2.7905, 2.5328, 3.2995, 2.7905, 2.8163, 4.2184, 3.8162, 2.4973, 3.2863, 3.0664, 3.8614,
			3.2146, 2.8163, 2.8163, 2.7931, 2.8335, 2.5328, 3.8162, 2.7905, 2.4973, 2.5328, 3.2863, 3.2995, 2.8936,
			3.0664, 2.7905, 3.0664, 2.8163, 3.8614, 2.8335, 2.7931, 2.8936, 2.4895, 3.2863, 2.8163, 2.4973, 2.8163,
			2.7931, 2.7905, 3.0597, 2.7931, 3.0627, 2.8936, 2.4973, 3.2995, 2.4973, 2.4895, 2.7931, 3.0664, 3.2995,
			2.5328, 3.2146, 2.7931, 2.7905, 3.2863, 3.2863, 2.8335, 2.4973, 2.8335, 3.0664, 2.8335, 3.2146, 2.8163,
			3.0597, 3.2995, 3.2146, 2.8283, 2.6607, 3.0627, 3.0627, 3.8162, 2.8163, 3.2995, 3.0664, 3.8614, 2.8335,
			3.0597, 3.8614, 3.2863, 2.5328, 2.7931, 2.7905, 3.0664, 3.2146, 2.4973
		};

		EnergyMap emap;
		for ( int ii = 1; ii <= 113; ++ii ) {
			emap.zero();
			paa_energy->residue_energy( the_pose->residue(ii), *the_pose, emap );
			TS_ASSERT_DELTA( emap[ p_aa ], correct_answers[ ii - 1 ], TOLERATED_ERROR );
		}

	}

};


