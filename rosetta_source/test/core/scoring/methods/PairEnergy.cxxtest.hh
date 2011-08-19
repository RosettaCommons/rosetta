// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/PairEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::PairEnergy.cc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

// Unit headers
#include <core/scoring/methods/PairEnergy.hh>

#include <platform/types.hh>

// Package Headers
#include <test/util/pose_funcs.hh>
#include <test/util/deriv_funcs.hh>
#include <test/core/init_util.hh>

#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <numeric/conversions.hh>

//Auto Headers
#include <core/conformation/Atom.fwd.hh>
#include <core/id/NamedStubID.hh>
#include <core/io/pdb/file_data.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/scoring/types.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <ObjexxFCL/FArray.fwd.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class PairEnergyTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void test_pairE_start_score_start_func_match_w_total_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_pair, 0.5 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.validate_start_func_matches_start_score( -0.524698508539912, false, 1e-6 );
	}

	void test_pairE_deriv_check_w_total_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_pair, 0.5 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( false, 1e-6 );
	}

	void test_pairE_start_score_start_func_match_w_partial_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_pair, 0.5 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.validate_start_func_matches_start_score( -0.524698508539912, false, 1e-6 );

	}

	void test_pairE_deriv_check_w_partial_flexibility()
	{
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_pair, 0.5 );
		kinematics::MoveMap movemap( create_trpcage_movemap_to_allow_bb10_freedom() );
		AtomDerivValidator adv;
		adv.set_pose( pose );
		adv.set_score_function( sfxn );
		adv.set_movemap( movemap );
		adv.simple_deriv_check( false, 1e-6 );
	}

};


