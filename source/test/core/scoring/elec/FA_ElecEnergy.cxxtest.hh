// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/elec/FA_ElecEnergy.cxxtest.hh
/// @brief  Test the score and derivative evaluation for the FA_ElecEnergy class
/// @author Andrew Leaver-Fay

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/util/pdb1ubq.hh>
#include <test/core/init_util.hh>

// Package headers
#include <core/scoring/elec/FA_ElecEnergy.hh>
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

class FA_ElecEnergyTests : public CxxTest::TestSuite {

public:
  void setUp() {
		core_init();
	}

	void tearDown(){}

	void test_elec_pre_talaris_2013_settings()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core_init_with_additional_options( "-restore_pre_talaris_2013_behavior -override_rsd_type_limit" );

		core::pose::Pose pose = create_twores_1ubq_pose();

		core::Real const TOL(1e-5);

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_elec, 1 );

		methods::EnergyMethodOptions options; // default is what we want

		core::scoring::etable::coulomb::Coulomb coulomb( options );
		//options.show(TR);

		EnergyMap emap;

		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												 pose.residue(2).xyz(pose.residue(1).atom_index("CB") ), -0.18 ), 0, TOL ); // CB-CB sits at 5.56174
		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												 pose.residue(2).xyz(pose.residue(1).atom_index("CA") ), 0.07 ), -0.006643886, TOL ); // CB-CA at 4.49784 (no count pair)
		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
																												 pose.residue(2).xyz(pose.residue(1).atom_index("N") ), -0.47 ), -3.175849739, TOL ); // C-N at 1.29914 (no count pair)

		core::scoring::elec::FA_ElecEnergy elec( options );
		elec.residue_pair_energy( pose.residue(1), pose.residue(2), pose, sfxn, emap );
		// Value calculated for 1.5/5.5 min/max, die=10r,
		// with no fa_elec on atoms 1,2,or 3 bonds apart, and a scaling of 0.2 on those 4 bonds apart
		TS_ASSERT_DELTA( emap[ fa_elec ] , 0.337019896208238, TOL);
		TS_ASSERT_DELTA( sfxn(pose), 0.337019896208238, TOL);

	}

	void test_elec_default_settings()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::Pose pose = create_twores_1ubq_pose();

		core::Real const TOL(1e-5);

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( fa_elec, 1 );

		methods::EnergyMethodOptions options; // default is what we want

		core::scoring::etable::coulomb::Coulomb coulomb( options );
		//options.show(TR);

		EnergyMap emap;

		//std::cout.precision( 6 );
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("CB") ), -0.18 ) <<
		//	", TOL ); // CB-CB sits at 5.56174" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"CA\") ), 0.07 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("CA") ), 0.07 ) <<
		//	", TOL ); // CB-CA at 4.49784 (no count pair)" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"C\") ), 0.51," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"N\") ), -0.47 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("N") ), -0.47 ) <<
		//	", TOL ); // C-N at 1.29914 (no count pair)" << std::endl;

		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												pose.residue(2).xyz(pose.residue(1).atom_index("CB") ), -0.18 ), 0, TOL ); // CB-CB sits at 5.56174
		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												pose.residue(2).xyz(pose.residue(1).atom_index("CA") ), 0.07 ), -0.00664392, TOL ); // CB-CA at 4.49784 (no count pair)
		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
																												pose.residue(2).xyz(pose.residue(1).atom_index("N") ), -0.47 ), -2.76037, TOL ); // C-N at 1.29914 (no count pair)

		core::scoring::elec::FA_ElecEnergy elec( options );
		elec.residue_pair_energy( pose.residue(1), pose.residue(2), pose, sfxn, emap );

		//std::cout.precision(16);
		//std::cout << "TS_ASSERT_DELTA( emap[ fa_elec ], " << emap[ fa_elec ] << ", TOL);" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( sfxn(pose)," <<  sfxn(pose) << ", TOL);" << std::endl;

		TS_ASSERT_DELTA( emap[ fa_elec ], 0.3312078368369487, TOL);
		TS_ASSERT_DELTA( sfxn(pose),0.3312078368369487, TOL);

	}


	void test_elec_max_dis()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::Pose pose = create_twores_1ubq_pose();

		core::Real const TOL(1e-5);

		methods::EnergyMethodOptions options;
		options.elec_max_dis( 7.0 );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_elec, 1 );

		core::scoring::etable::coulomb::Coulomb coulomb( options );
		options.show(TR);

		EnergyMap emap;

		//std::cout.precision( 6 );
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("CB") ), -0.18 ) <<
		//	", TOL ); // CB-CB sits at 5.56174" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"CA\") ), 0.07 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("CA") ), 0.07 ) <<
		//	", TOL ); // CB-CA at 4.49784 (no count pair)" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"C\") ), 0.51," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"N\") ), -0.47 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("N") ), -0.47 ) <<
		//	", TOL ); // C-N at 1.29914 (no count pair)" << std::endl;

		//TR << "Distance : " << pose.residue(1).xyz(pose.residue(1).atom_index("CB")).distance(pose.residue(2).xyz(pose.residue(2).atom_index("CB"))) << std::endl;

		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												pose.residue(2).xyz(pose.residue(1).atom_index("CB") ), -0.18 ), 0.0115767, TOL ); // CB-CB sits at 5.56174
		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												pose.residue(2).xyz(pose.residue(1).atom_index("CA") ), 0.07 ), -0.0117772, TOL ); // CB-CA at 4.49784 (no count pair)
		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
																												pose.residue(2).xyz(pose.residue(1).atom_index("N") ), -0.47 ), -2.85802, TOL ); // C-N at 1.29914 (no count pair)

		core::scoring::elec::FA_ElecEnergy elec( options );

		elec.residue_pair_energy( pose.residue(1), pose.residue(2), pose, sfxn, emap );
		// Value calculated for 1.5/7.0 min/max, die=10r,
		// with no fa_elec on atoms 1,2,or 3 bonds apart, and a scaling of 0.2 on those 4 bonds apart

		//std::cout.precision(16);
		//std::cout << "TS_ASSERT_DELTA( emap[ fa_elec ], " << emap[ fa_elec ] << ", TOL);" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( sfxn(pose)," <<  sfxn(pose) << ", TOL);" << std::endl;

		TS_ASSERT_DELTA( emap[ fa_elec ], 0.62052938253723, TOL);
		TS_ASSERT_DELTA( sfxn(pose),0.62052938253723, TOL);

	}

	void test_elec_min_dis()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::Pose pose = create_twores_1ubq_pose();

		core::Real const TOL(1e-5);

		methods::EnergyMethodOptions options;
		options.elec_min_dis( 3.0 );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_elec, 1 );

		core::scoring::etable::coulomb::Coulomb coulomb( options );
		//options.show(TR);

		EnergyMap emap;

		//std::cout.precision( 6 );
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("CB") ), -0.18 ) <<
		//	", TOL ); // CB-CB sits at 5.56174" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"CA\") ), 0.07 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("CA") ), 0.07 ) <<
		//	", TOL ); // CB-CA at 4.49784 (no count pair)" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"C\") ), 0.51," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"N\") ), -0.47 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("N") ), -0.47 ) <<
		//	", TOL ); // C-N at 1.29914 (no count pair)" << std::endl;

		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												pose.residue(2).xyz(pose.residue(1).atom_index("CB") ), -0.18 ), 0, TOL ); // CB-CB sits at 5.56174
		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												pose.residue(2).xyz(pose.residue(1).atom_index("CA") ), 0.07 ), -0.00664392, TOL ); // CB-CA at 4.49784 (no count pair)
		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
																												pose.residue(2).xyz(pose.residue(1).atom_index("N") ), -0.47 ), -0.602561, TOL ); // C-N at 1.29914 (no count pair)

		core::scoring::elec::FA_ElecEnergy elec( options );
		elec.residue_pair_energy( pose.residue(1), pose.residue(2), pose, sfxn, emap );

		// Value calculated for 3.0/5.5 min/max, die=10r,
		// with no fa_elec on atoms 1,2,or 3 bonds apart, and a scaling of 0.2 on those 4 bonds apart
		//std::cout.precision(16);
		//std::cout << "TS_ASSERT_DELTA( emap[ fa_elec ], " << emap[ fa_elec ] << ", TOL);" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( sfxn(pose)," <<  sfxn(pose) << ", TOL);" << std::endl;

		TS_ASSERT_DELTA( emap[ fa_elec ], 0.3191849378777624, TOL);
		TS_ASSERT_DELTA( sfxn(pose),0.3191849378777624, TOL);

	}

	void test_elec_die()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::Pose pose = create_twores_1ubq_pose();

		core::Real const TOL(1e-5);

		methods::EnergyMethodOptions options;
		options.elec_die( 4.0 );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_elec, 1 );

		core::scoring::etable::coulomb::Coulomb coulomb( options );
		//options.show(TR);

		EnergyMap emap;

		//std::cout.precision( 6 );
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("CB") ), -0.18 ) <<
		//	", TOL ); // CB-CB sits at 5.56174" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"CA\") ), 0.07 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("CA") ), 0.07 ) <<
		//	", TOL ); // CB-CA at 4.49784 (no count pair)" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"C\") ), 0.51," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"N\") ), -0.47 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("N") ), -0.47 ) <<
		//	", TOL ); // C-N at 1.29914 (no count pair)" << std::endl;


		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												pose.residue(2).xyz(pose.residue(1).atom_index("CB") ), -0.18 ), 0, TOL ); // CB-CB sits at 5.56174
		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												pose.residue(2).xyz(pose.residue(1).atom_index("CA") ), 0.07 ), -0.0166098, TOL ); // CB-CA at 4.49784 (no count pair)
		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
																												pose.residue(2).xyz(pose.residue(1).atom_index("N") ), -0.47 ), -6.90093, TOL ); // C-N at 1.29914 (no count pair)

		core::scoring::elec::FA_ElecEnergy elec( options );
		elec.residue_pair_energy( pose.residue(1), pose.residue(2), pose, sfxn, emap );
		// Value calculated for 1.5/5.5 min/max, die=4r,
		// with no fa_elec on atoms 1,2,or 3 bonds apart, and a scaling of 0.2 on those 4 bonds apart

		//std::cout.precision(16);
		//std::cout << "TS_ASSERT_DELTA( emap[ fa_elec ], " << emap[ fa_elec ] << ", TOL);" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( sfxn(pose)," <<  sfxn(pose) << ", TOL);" << std::endl;

		TS_ASSERT_DELTA( emap[ fa_elec ], 0.8280195920923709, TOL);
		TS_ASSERT_DELTA( sfxn(pose),0.8280195920923709, TOL);

	}

	void test_elec_no_dis_dep_die()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		core::pose::Pose pose = create_twores_1ubq_pose();

		core::Real const TOL(1e-5);

		methods::EnergyMethodOptions options;
		options.elec_no_dis_dep_die( true );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_elec, 1 );

		core::scoring::etable::coulomb::Coulomb coulomb( options );
		//options.show(TR);

		EnergyMap emap;

		//std::cout.precision( 6 );
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("CB") ), -0.18 ) <<
		//	", TOL ); // CB-CB sits at 5.56174" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"CB\") ), -0.18," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"CA\") ), 0.07 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("CA") ), 0.07 ) <<
		//	", TOL ); // CB-CA at 4.49784 (no count pair)" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index(\"C\") ), 0.51," << std::endl;
		//std::cout << "pose.residue(2).xyz(pose.residue(1).atom_index(\"N\") ), -0.47 ), " <<
		//	coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
		//	pose.residue(2).xyz(pose.residue(1).atom_index("N") ), -0.47 ) <<
		//	", TOL ); // C-N at 1.29914 (no count pair)" << std::endl;


		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												pose.residue(2).xyz(pose.residue(1).atom_index("CB") ), -0.18 ), 0, TOL ); // CB-CB sits at 5.56174
		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("CB") ), -0.18,
																												pose.residue(2).xyz(pose.residue(1).atom_index("CA") ), 0.07 ), -0.0164394, TOL ); // CB-CA at 4.49784 (no count pair)
		TS_ASSERT_DELTA( coulomb.eval_atom_atom_fa_elecE( pose.residue(1).xyz(pose.residue(1).atom_index("C") ), 0.51,
																												pose.residue(2).xyz(pose.residue(1).atom_index("N") ), -0.47 ), -3.4213, TOL ); // C-N at 1.29914 (no count pair)


		core::scoring::elec::FA_ElecEnergy elec( options );
		elec.residue_pair_energy( pose.residue(1), pose.residue(2), pose, sfxn, emap );
		// Value calculated for 1.5/5.5 min/max, die=10 (no r dependance),
		// with no fa_elec on atoms 1,2,or 3 bonds apart, and a scaling of 0.2 on those 4 bonds apart

		//std::cout.precision(16);
		//std::cout << "TS_ASSERT_DELTA( emap[ fa_elec ], " << emap[ fa_elec ] << ", TOL);" << std::endl;
		//std::cout << "TS_ASSERT_DELTA( sfxn(pose)," <<  sfxn(pose) << ", TOL);" << std::endl;

		TS_ASSERT_DELTA( emap[ fa_elec ], 0.834303865499549, TOL);
		TS_ASSERT_DELTA( sfxn(pose),0.834303865499549, TOL);

	}

	/// @brief Setup for minimization using a move map that says "minimize all bb and sc torsions"
	/// Make sure that start_score matches start_func.
	void test_elec_deriv_check_w_full_torsional_flexibility()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = pdb1ubq5to13_pose();
		ScoreFunction sfxn;
		sfxn.set_weight( fa_elec, 0.75 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );

		/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
		/// that the analytic norm matches the numeric norm to 1e-3.
		adv.simple_deriv_check( true, 1e-6 );

	}

	/// @brief Smoothed fa-elec derivative check.
	/// Setup for minimization using a move map that says "minimize all bb and sc torsions"
	/// Make sure that start_score matches start_func.
	void test_elec_deriv_check_w_full_torsional_flexibility_smoothed()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = pdb1ubq5to13_pose();
		methods::EnergyMethodOptions options;
		options.smooth_fa_elec( true );

		ScoreFunction sfxn;
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_elec, 0.75 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );

		/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
		/// that the analytic norm matches the numeric norm to 1e-3.
		adv.simple_deriv_check( true, 1e-6 );

	}

	/// @brief Smoothed fa-elec derivative check.
	/// Setup for minimization using a move map that says "minimize all bb and sc torsions"
	/// Make sure that start_score matches start_func.
	void test_elec_deriv_check_w_full_torsional_flexibility_smoothed_wo_ddd()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = pdb1ubq5to13_pose();
		methods::EnergyMethodOptions options;
		options.smooth_fa_elec( true );
		options.elec_no_dis_dep_die( true );

		ScoreFunction sfxn;
		sfxn.set_energy_method_options( options );
		sfxn.set_weight( fa_elec, 0.75 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );

		/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
		/// that the analytic norm matches the numeric norm to 1e-3.
		adv.simple_deriv_check( true, 1e-6 );

	}


	/// @brief Smoothed fa-elec finalize energy check.
	void test_elec_nblist_autoupdate_no_smoothing()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = pdb1ubq5to13_pose();
		methods::EnergyMethodOptions em_options;
		em_options.smooth_fa_elec( false );

		ScoreFunction sfxn;
		sfxn.set_energy_method_options( em_options );
		sfxn.set_weight( fa_elec, 0.75 );

		kinematics::MoveMap movemap;
		movemap.set_bb( true );
		movemap.set_chi( true );

		AtomDerivValidator adv( pose, sfxn, movemap );
		adv.set_nblist_auto_update( true );

		/// This call runs a numeric deriv check on all the free dofs in the system and makes sure
		/// that the analytic norm matches the numeric norm to 1e-3.
		adv.simple_deriv_check( true, 1e-6 );

	}

	/// @brief Smoothed fa-elec finalize energy check.
	void test_elec_nblist_autoupdate_with_smoothing()
	{
		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::optimization;

		Pose pose = pdb1ubq5to13_pose();
		methods::EnergyMethodOptions em_options;
		em_options.smooth_fa_elec( true );

		ScoreFunction sfxn;
		sfxn.set_energy_method_options( em_options );
		sfxn.set_weight( fa_elec, 0.75 );

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
