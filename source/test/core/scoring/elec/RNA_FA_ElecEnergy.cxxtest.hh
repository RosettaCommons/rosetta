// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/elec/RNA_FA_ElecEnergy.cxxtest.hh
/// @brief  Test the score and derivative evaluation for the RNA_FA_ElecEnergy class
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1lnt.hh>
#include <test/core/init_util.hh>

// Package headers
#include <core/scoring/elec/RNA_FA_ElecEnergy.hh>
#include <core/scoring/etable/coulomb/Coulomb.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyMap.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("core.scoring.elec.FA_ElecEnergy.cxxtest");

using namespace core;

class RNA_FA_ElecEnergyTests : public CxxTest::TestSuite {

public:
	void setUp() {
		core_init();
	}

	void tearDown(){}

	void test_score_rna_pdb() {
		core::pose::Pose rna_pose = pdb1lnt_pose();
		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::fa_elec_rna_phos_phos, 0.75 );
		core::Real const score = sfxn( rna_pose );
		//std::cout.precision(16);
		//std::cout << "score rna pdb: " << score << std::endl;
		TS_ASSERT_DELTA( score, -0.3867971789450972, 1e-6 );
	}

	/// @brief Smoothed fa-elec derivative check.
	/// Setup for minimization using a move map that says "minimize all bb and sc torsions"
	/// Make sure that start_score matches start_func.
	void test_rnaelec_deriv_check_w_full_torsional_flexibility_smoothed_wo_ddd()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = pdb1lnt_pose();
		methods::EnergyMethodOptions options;
		options.smooth_fa_elec( true );
		options.elec_no_dis_dep_die( true );

		ScoreFunction sfxn;
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_elec_rna_phos_phos, 0.75 );
		/*sfxn.set_weight( fa_elec_rna_phos_base, 0.75 );
		sfxn.set_weight( fa_elec_rna_phos_sugr, 0.75 );
		sfxn.set_weight( fa_elec_rna_base_base, 0.75 );
		sfxn.set_weight( fa_elec_rna_sugr_base, 0.75 );
		sfxn.set_weight( fa_elec_rna_sugr_sugr, 0.75 );
		*/
		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );

		/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
		/// that the analytic norm matches the numeric norm to 1e-3.
		adv.simple_deriv_check( true, 1e-6 );

	}


	/// @brief Smoothed fa-elec finalize energy check.
	void test_rnaelec_nblist_autoupdate_no_smoothing()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = pdb1lnt_pose();
		methods::EnergyMethodOptions em_options;
		em_options.smooth_fa_elec( false );

		ScoreFunction sfxn;
		sfxn.set_energy_method_options( em_options );
		sfxn.set_weight( fa_elec_rna_phos_phos, 0.75 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.set_nblist_auto_update( true );

		/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
		/// that the analytic norm matches the numeric norm to 1e-3.
		adv.simple_deriv_check( true, 1e-6 );

	}

};
