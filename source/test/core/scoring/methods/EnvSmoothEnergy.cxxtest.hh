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
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <platform/types.hh>

// Unit headers
#include <core/scoring/methods/EnvSmoothEnergy.hh>

// Numeric headers
#include <numeric/conversions.hh>


// Project headers
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

//Auto Headers
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class EnvSmoothEnergyTests : public CxxTest::TestSuite {

public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		using namespace core;
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_setup_for_minimizing_with_envsmooth_partial_bb_flex_and_full_sc_flex()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn; sfxn.set_weight( envsmooth,   0.4 );

		kinematics::MoveMap movemap; movemap.set_bb( 10, true ); movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.validate_start_func_matches_start_score();
	}

	void dont_test_atom_tree_minimize_with_envsmooth_and_etable()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_atr,    0.5 );
		sfxn.set_weight( fa_rep,    0.5 );
		sfxn.set_weight( envsmooth, 0.4 );

		AtomTreeMinimizer minimizer;
		Real start_score = sfxn(pose);
		//std::cout << "start score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -27.89673840790607, start_score, 1e-12 );


		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );

		minimizer.run( pose, movemap, sfxn, min_options );

		Real end_score = sfxn(pose);
		//std::cout << "end score: " << sfxn(pose) << std::endl;
		TS_ASSERT_DELTA( -28.10526841844615, end_score, 1e-12 );
	}

	/// @brief Create a move map that has the central residue changing, while everything
	/// else is held fixed; then the domain map will have color 1 for residues 1-9, and
	/// color 2 for residues 11-20.
	void test_setup_for_minimizing_with_envsmooth_and_partial_bbflex()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn; sfxn.set_weight( envsmooth,   0.4 );

		kinematics::MoveMap movemap; movemap.set_bb( 10, true );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.validate_start_func_matches_start_score();
	}

	void test_envsmooth_deriv_check_w_full_bb_flexibility()
	{
		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn; sfxn.set_weight( envsmooth, 0.4 );
		kinematics::MoveMap movemap; movemap.set_bb( true );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6 );
	}


};


