// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/RamachandranEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::RamachandranEnergy.cc
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/methods/YHHPlanarityEnergy.hh>

#include <platform/types.hh>

// Package Headers
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>
#include <test/core/init_util.hh>

// AUTO-REMOVED #include <core/id/DOF_ID.hh>
// AUTO-REMOVED #include <core/id/TorsionID.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <basic/options/option.hh>

// AUTO-REMOVED #include <core/optimization/MinimizerOptions.hh>
// AUTO-REMOVED #include <core/optimization/AtomTreeMinimizer.hh>

#include <numeric/conversions.hh>

//Auto Headers
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class YHHPlanarityEnergyTests : public CxxTest::TestSuite {

public:

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
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_yhh_planarity_start_score_start_func_match_w_total_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( yhh_planarity, 0.5 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( 0.06583992062591906, false /* print start score?*/ );

	}


	void test_yhh_planarity_start_score_start_func_match_w_partial_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( yhh_planarity, 0.5 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( 0.06583992062591906,  false /* print start score?*/ );

	}

	/// Ensure norm matches norm-numeric
	void test_yhh_planarity_deriv_check_w_total_flexibility()
	{
		using namespace core::kinematics;
		using namespace core::pose;
		using namespace core::scoring;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( yhh_planarity, 0.5 );
		MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6 );
	}

};


