// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constaints/AngleConstraint.cxxtest.hh
/// @brief  test suite for angle constraints
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>

// AUTO-REMOVED #include <core/conformation/Residue.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

// Project headers
// Auto-header: duplicate removed #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerMap.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerOptions.hh>
// AUTO-REMOVED #include <core/optimization/AtomTreeMinimizer.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

// AUTO-REMOVED #include <numeric/conversions.hh>
#include <numeric/constants.hh>

//Auto Headers
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.scoring.constraints.DihedralConstraint.cxxtest");

using namespace core;

class ConstraintsEnergyTests : public CxxTest::TestSuite
{

public:
	ConstraintsEnergyTests() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_atom_tree_minimize_with_dist_cst()
	{
		using namespace core;
		using namespace core::id;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( atom_pair_constraint, 0.5 );

		// Try to draw the tyrosine 3 sidechain towards the backbone oxygen on residue 19.
		FuncOP harmonic_func = new HarmonicFunc( 4.0, 1.0 );
		AtomPairConstraintOP cst1 = new AtomPairConstraint(
			AtomID( pose.residue( 3  ).atom_index( "CE2" ), 3 ),
			AtomID( pose.residue( 19 ).atom_index( "O"   ), 19),
			harmonic_func );
		pose.add_constraint( cst1 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6 );
	}

	void test_atom_tree_minimize_with_dist_cst_in_fixed_domain()
	{
		using namespace core;
		using namespace core::id;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( atom_pair_constraint, 0.5 );

		// Try to push the tyrosine 3 sidechain away from CZ3.
		FuncOP harmonic_func = new HarmonicFunc( 5.0, 1.0 );
		AtomPairConstraintOP cst1 = new AtomPairConstraint(
			AtomID( pose.residue( 3  ).atom_index( "CZ" ), 3 ),
			AtomID( pose.residue( 6 ).atom_index(  "CZ3"   ), 6),
			harmonic_func );
		pose.add_constraint( cst1 );

		kinematics::MoveMap movemap;
		movemap.set_bb( 10, true );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-6 );
	}

	void test_atom_tree_minimize_with_intraresidue_dist_cst()
	{
		using namespace core;
		using namespace core::id;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( atom_pair_constraint, 0.5 );

		// Try to draw the tyrosine 3 sidechain towards the backbone oxygen on residue 19.
		// of course, this pair is a fixed distance in this minimization routine given the domain map
		FuncOP harmonic_func = new HarmonicFunc( 3.5, 1.0 );
		AtomPairConstraintOP cst1 = new AtomPairConstraint(
			AtomID( pose.residue( 2 ).atom_index( "CD2" ), 2 ),
			AtomID( pose.residue( 2 ).atom_index( "N"   ), 2 ),
			harmonic_func );
		pose.add_constraint( cst1 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-5 );

	}

	void test_atom_tree_minimize_with_interresidue_angle_cst()
	{
		using namespace core;
		using namespace core::id;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::optimization;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( angle_constraint, 0.5 );

		FuncOP harmonic_func = new HarmonicFunc( numeric::constants::d::pi_over_2, numeric::constants::d::deg2rad * 5 );
		ConstraintOP cst1 = new AngleConstraint(
			AtomID( pose.residue( 6  ).atom_index( "CZ3" ), 6 ),
			AtomID( pose.residue( 3 ).atom_index(  "CZ"  ), 3 ),
			AtomID( pose.residue( 3 ).atom_index(  "CG"  ), 3 ),
			harmonic_func );
		pose.add_constraint( cst1 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-5 );
	}

	void test_atom_tree_minimize_with_threebody_angle_cst()
	{
		using namespace core;
		using namespace core::id;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::constraints;

		Pose pose = create_trpcage_ideal_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( angle_constraint, 0.5 );

		FuncOP harmonic_func = new HarmonicFunc( numeric::constants::d::pi_over_2, numeric::constants::d::deg2rad * 5 );
		ConstraintOP cst1 = new AngleConstraint(
			AtomID( pose.residue( 6  ).atom_index( "CZ3" ),  6 ),
			AtomID( pose.residue( 3  ).atom_index( "CZ"  ),  3 ),
			AtomID( pose.residue( 19 ).atom_index( "O"   ), 19 ),
			harmonic_func );
		pose.add_constraint( cst1 );
		kinematics::MoveMap movemap( create_movemap_to_allow_all_torsions() );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.simple_deriv_check( true, 1e-5 );

	}
};
