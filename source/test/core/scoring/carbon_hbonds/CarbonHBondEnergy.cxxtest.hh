// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/CarbonHBondEnergy.cxxtest.hh
/// @brief  Test the score and derivative evaluation for the CarbonHBondEnergy class
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1ubq.hh>
#include <test/core/init_util.hh>

// Package headers
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerMap.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerOptions.hh>
// AUTO-REMOVED #include <core/optimization/AtomTreeMinimizer.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>

// AUTO-REMOVED #include <basic/Tracer.hh>

// AUTO-REMOVED #include <core/id/DOF_ID.hh>
// AUTO-REMOVED #include <core/id/TorsionID.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <basic/options/option.hh>


// Utility headers
#include <utility/vector1.hh>

//Auto Headers


using namespace core;

class CarbonHBondEnergyTests : public CxxTest::TestSuite {

public:
  void setUp() {
		core_init_with_additional_options( "-restore_pre_talaris_2013_behavior -override_rsd_type_limit");
	}

	void tearDown(){}

	/// @brief Setup for minimization using a move map that says "minimize all bb and sc torsions"
	/// Make sure that start_score matches start_func.
	void test_chbonds_deriv_check_w_full_torsional_flexibility()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = pdb1ubq5to13_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( ch_bond_bb_bb, 0.75 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );

		/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
		/// that the analytic norm matches the numeric norm to 1e-3.
		adv.simple_deriv_check( true, 2e-1 );

	}

	/// @brief The deriv check above is very low resolution (i.e. wrong) but the derivatives
	/// I'm now adding match the derivatives observed in trunk.  When the derivatives
	/// are fixed, remove this test.
	void test_chbonds_deriv_check_matches_trunk_r38117()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = pdb1ubq5to13_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( ch_bond_bb_bb, 0.75 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.set_nonzero_deriv_only( true );

		//adv.compute_pose_atom_derivs(); <--- this function call generated the f1/f2s below

		using namespace core;
		using namespace core::id;
		AtomDerivList adl;
		adl.add( AtomDeriv( AtomID( 4, 1), 23.04123023485344,-8.495397081825313,-30.31768491462694,-0.03022263897114251,0.9899924456287443,-0.3003773455033463));
		adl.add( AtomDeriv( AtomID( 11, 2), 7.197844267884681,-0.6801661037170241,-14.5096977067363,-0.7712133151326138,-0.5398437833845259,-0.3572707027849681));
		adl.add( AtomDeriv( AtomID( 4, 7), -6.334790502302728,-0.2540206135385876,14.05826504628001,0.7712133151326138,0.5398437833845259,0.3572707027849681));
		adl.add( AtomDeriv( AtomID( 9, 8), -23.04123023485344,8.495397081825313,30.31768491462695,0.03022263897114251,-0.9899924456287443,0.3003773455033463));

		adv.validate_atom_deriv_list( adl );

	}


};
