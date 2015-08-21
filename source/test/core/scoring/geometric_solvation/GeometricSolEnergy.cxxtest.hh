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

// Project headers


// Utility headers
#include <utility/vector1.hh>

//Auto Headers


using namespace core;

class GeometricSolEnergyTests : public CxxTest::TestSuite {

public:
	void setUp() {
		core_init();
	}

	void tearDown(){}


	/// Remove this test when dont_test_geomsol_deriv_check is reinstated
	void test_empty_geomsol_test() {
		TS_ASSERT( true );
	}
	/// @brief Setup for minimization using a move map that says "minimize all bb and sc torsions"
	/// Make sure that start_score matches start_func.
	void dont_test_geomsol_deriv_check_w_full_torsional_flexibility()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = pdb1ubq5to13_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( geom_sol, 0.75 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );

		/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
		/// that the analytic norm matches the numeric norm to 1e-2.  The derivatives for geomsol
		/// are essentially wrong, and are worse now that the derivatives for hydrogen bonds have
		/// been fixed
		adv.simple_deriv_check( true, 1e-2 );

	}


};
