// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hackelec/HackElecEnergy.cxxtest.hh
/// @brief  Test the score and derivative evaluation for the HackElecEnergy class
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1ubq.hh>
#include <test/core/init_util.hh>

// Package headers
#include <core/scoring/hackelec/HackElecEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("core.scoring.hackelec.HackElecEnergy.cxxtest");

using namespace core;

class HackElecEnergyTests : public CxxTest::TestSuite {

public:
  void setUp() {
		core_init();
	}

	void tearDown(){}

	void test_hackelec_default_settings()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::Pose pose = create_twores_1ubq_pose();

		core::Real const TOL(1e-5);

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( hack_elec, 1 );

		methods::EnergyMethodOptions options; // default is what we want

		core::scoring::hackelec::HackElecEnergy hackelec( options );
		//options.show(TR);

		EnergyMap emap;

		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												 pose.residue(2).xyz(pose.residue(1).atom_index("CB") ), -0.18 ), 0, TOL ); // CB-CB sits at 5.56174
		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												 pose.residue(2).xyz(pose.residue(1).atom_index("CA") ), 0.07 ), -0.006643886, TOL ); // CB-CA at 4.49784 (no count pair)
		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
																												 pose.residue(2).xyz(pose.residue(1).atom_index("N") ), -0.47 ), -3.175849739, TOL ); // C-N at 1.29914 (no count pair)

		hackelec.residue_pair_energy( pose.residue(1), pose.residue(2), pose, sfxn, emap );
		// Value calculated for 1.5/5.5 min/max, die=10r,
		// with no hack_elec on atoms 1,2,or 3 bonds apart, and a scaling of 0.2 on those 4 bonds apart
		TS_ASSERT_DELTA( emap[ hack_elec ] , 0.337019896208238, TOL);
		TS_ASSERT_DELTA( sfxn(pose), 0.337019896208238, TOL);

	}

	void test_hackelec_max_dis()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::Pose pose = create_twores_1ubq_pose();

		core::Real const TOL(1e-5);

		methods::EnergyMethodOptions options;
		options.hackelec_max_dis( 7.0 );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( hack_elec, 1 );

		core::scoring::hackelec::HackElecEnergy hackelec( options );
		options.show(TR);

		EnergyMap emap;

		TR << "Distance : " << pose.residue(1).xyz(pose.residue(1).atom_index("CB")).distance(pose.residue(2).xyz(pose.residue(2).atom_index("CB"))) << std::endl;
		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												 pose.residue(2).xyz(pose.residue(2).atom_index("CB") ), -0.18 ), 0.012438172, TOL ); // CB-CB sits at 5.56174
		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												 pose.residue(2).xyz(pose.residue(2).atom_index("CA") ), 0.07 ), -0.011777133, TOL ); // CB-CA at 4.49784 (no count pair)
		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
																												 pose.residue(2).xyz(pose.residue(2).atom_index("N") ), -0.47 ), -3.273503647, TOL ); // C-N at 1.29914 (no count pair)

		hackelec.residue_pair_energy( pose.residue(1), pose.residue(2), pose, sfxn, emap );
		// Value calculated for 1.5/7.0 min/max, die=10r,
		// with no hack_elec on atoms 1,2,or 3 bonds apart, and a scaling of 0.2 on those 4 bonds apart
		TS_ASSERT_DELTA( emap[ hack_elec ] , 0.589919643, TOL);
		TS_ASSERT_DELTA( sfxn(pose), 0.589919643, TOL);
	}

	void test_hackelec_min_dis()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::Pose pose = create_twores_1ubq_pose();

		core::Real const TOL(1e-5);

		methods::EnergyMethodOptions options;
		options.hackelec_min_dis( 3.0 );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( hack_elec, 1 );

		core::scoring::hackelec::HackElecEnergy hackelec( options );
		//options.show(TR);

		EnergyMap emap;

		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												 pose.residue(2).xyz(pose.residue(2).atom_index("CB") ), -0.18 ), 0, TOL ); // CB-CB sits at 5.56174
		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												 pose.residue(2).xyz(pose.residue(2).atom_index("CA") ), 0.07 ), -0.006643886, TOL ); // CB-CA at 4.49784 (no count pair)
		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
																												 pose.residue(2).xyz(pose.residue(2).atom_index("N") ), -0.47 ), -0.602560776, TOL ); // C-N at 1.29914 (no count pair)

		hackelec.residue_pair_energy( pose.residue(1), pose.residue(2), pose, sfxn, emap );
		// Value calculated for 3.0/5.5 min/max, die=10r,
		// with no hack_elec on atoms 1,2,or 3 bonds apart, and a scaling of 0.2 on those 4 bonds apart
		TS_ASSERT_DELTA( emap[ hack_elec ] , 0.319910941, TOL);
		TS_ASSERT_DELTA( sfxn(pose), 0.319910941, TOL);
	}

	void test_hackelec_die()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::Pose pose = create_twores_1ubq_pose();

		core::Real const TOL(1e-5);

		methods::EnergyMethodOptions options;
		options.hackelec_die( 4.0 );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( hack_elec, 1 );

		core::scoring::hackelec::HackElecEnergy hackelec( options );
		//options.show(TR);

		EnergyMap emap;

		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												 pose.residue(2).xyz(pose.residue(2).atom_index("CB") ), -0.18 ), 0, TOL ); // CB-CB sits at 5.56174
		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												 pose.residue(2).xyz(pose.residue(2).atom_index("CA") ), 0.07 ), -0.016609716, TOL ); // CB-CA at 4.49784 (no count pair)
		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
																												 pose.residue(2).xyz(pose.residue(2).atom_index("N") ), -0.47 ), -7.939624349, TOL ); // C-N at 1.29914 (no count pair)

		hackelec.residue_pair_energy( pose.residue(1), pose.residue(2), pose, sfxn, emap );
		// Value calculated for 1.5/5.5 min/max, die=4r,
		// with no hack_elec on atoms 1,2,or 3 bonds apart, and a scaling of 0.2 on those 4 bonds apart
		TS_ASSERT_DELTA( emap[ hack_elec ] , 0.842549741, TOL);
		TS_ASSERT_DELTA( sfxn(pose), 0.842549741, TOL);
	}

	void test_hackelec_no_dis_dep_die()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::Pose pose = create_twores_1ubq_pose();

		core::Real const TOL(1e-5);

		methods::EnergyMethodOptions options;
		options.hackelec_no_dis_dep_die( true );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( hack_elec, 1 );

		core::scoring::hackelec::HackElecEnergy hackelec( options );
		//options.show(TR);

		EnergyMap emap;

		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												 pose.residue(2).xyz(pose.residue(2).atom_index("CB") ), -0.18 ), 0, TOL ); // CB-CB sits at 5.56174
		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												 pose.residue(2).xyz(pose.residue(2).atom_index("CA") ), 0.07 ), -0.016439276, TOL ); // CB-CA at 4.49784 (no count pair)
		TS_ASSERT_DELTA( hackelec.eval_atom_atom_hack_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
																												 pose.residue(2).xyz(pose.residue(2).atom_index("N") ), -0.47 ), -3.742965764, TOL ); // C-N at 1.29914 (no count pair)

		hackelec.residue_pair_energy( pose.residue(1), pose.residue(2), pose, sfxn, emap );
		// Value calculated for 1.5/5.5 min/max, die=10 (no r dependance),
		// with no hack_elec on atoms 1,2,or 3 bonds apart, and a scaling of 0.2 on those 4 bonds apart
		TS_ASSERT_DELTA( emap[ hack_elec ] , 0.850238278, TOL);
		TS_ASSERT_DELTA( sfxn(pose), 0.850238278, TOL);
	}

	/// @brief Setup for minimization using a move map that says "minimize all bb and sc torsions"
	/// Make sure that start_score matches start_func.
	void test_hackelec_deriv_check_w_full_torsional_flexibility()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = pdb1ubq5to13_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( hack_elec, 0.75 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );

		/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
		/// that the analytic norm matches the numeric norm to 1e-3.
		adv.simple_deriv_check( true, 1e-6 );

	}

};
