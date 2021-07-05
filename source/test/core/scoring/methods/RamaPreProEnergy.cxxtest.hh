// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/energy_methods/RamaPreProEnergy.cxxtest.hh
/// @brief  test suite for core::energy_methods::RamaPreProEnergy.cc
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/deriv_funcs.hh>
#include <test/util/cart_deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>


// Unit headers
#include <core/energy_methods/RamaPreProEnergy.hh>


// Package Headers
#include <core/id/PartialAtomID.hh>




// Project headers

//Auto Headers
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;

class RamaPreProEnergyTests : public CxxTest::TestSuite {

public:

	PoseOP the_pose;
	//RamachandranEnergyOP rama_energy;


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

	void test_atoms_w_dof_deriv() {
		Pose pose = create_trpcage_ideal_pose();
		core::energy_methods::RamaPreProEnergy rppe;

		auto atoms = rppe.atoms_with_dof_derivatives( pose.residue(10), pose );
		TS_ASSERT_EQUALS( atoms.size(), 5 );
		TS_ASSERT( std::find( atoms.begin(), atoms.end(), id::PartialAtomID(1, 10)) != atoms.end() );
		TS_ASSERT( std::find( atoms.begin(), atoms.end(), id::PartialAtomID(2, 10)) != atoms.end() );
		TS_ASSERT( std::find( atoms.begin(), atoms.end(), id::PartialAtomID(3, 10)) != atoms.end() );
		TS_ASSERT( std::find( atoms.begin(), atoms.end(), id::PartialAtomID(2, 9, 0)) != atoms.end() );
		TS_ASSERT( std::find( atoms.begin(), atoms.end(), id::PartialAtomID(1, 11, 0)) != atoms.end() );
	}


	// --------------- Test Cases --------------- //
	void test_setup_for_minimizing_with_rama_and_full_bb_flex()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;


		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn; sfxn.set_weight( rama_prepro,   0.4 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.validate_start_func_matches_start_score();
	}


	/// @brief Create a move map that has the central residue changing, while everything
	/// else is held fixed; then the domain map will have color 1 for residues 1-9, and
	/// color 2 for residues 11-20.
	void test_setup_for_minimizing_with_rama_and_partial_bbflex()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn; sfxn.set_weight( rama_prepro, 0.4 );

		kinematics::MoveMap movemap; movemap.set_bb( 10, true );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.validate_start_func_matches_start_score();
	}


	void test_rama_deriv_check_full_flexibility()
	{
		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn; sfxn.set_weight( rama_prepro, 0.4 );
		kinematics::MoveMap movemap; movemap.set_bb( true );
		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6 );
	}

	void test_rama_deriv_check_partial_flexibility()
	{
		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn; sfxn.set_weight( rama_prepro, 0.4 );

		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );
		movemap.set_bb( 11, true );
		movemap.set_bb( 12, true );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6 );
	}


	void test_rama_cart_deriv_check_partial_flexibility()
	{
		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn; sfxn.set_weight( rama_prepro, 0.4 );
		sfxn(pose);

		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );
		movemap.set_bb( 11, true );
		movemap.set_bb( 13, true );

		CartAtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-5 );
	}

	void test_symm_rama_cart_deriv_check_partial_flexibility()
	{
		core_init_with_additional_options( "-symmetry:symmetry_definition core/scoring/symmetry/sym_def.dat" );
		Pose pose;
		core::import_pose::pose_from_file(
			pose, "core/scoring/symmetry/test_in.pdb", core::import_pose::PDB_file );
		core::pose::symmetry::make_symmetric_pose( pose );

		ScoreFunction sfxn; sfxn.set_weight( rama_prepro, 0.4 );

		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );
		movemap.set_bb( 11, true );
		movemap.set_bb( 12, true );
		movemap.set_bb( 14, true );

		CartAtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-5 );
	}

};


