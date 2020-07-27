// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
//
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// /// @file test/core/scoring/HRFDynamics.cxxtest.hh
// /// @brief Unit tests for HRFDynamics score term.
// /// @details This is a unit test for score function HRFDynamics to ensure it's working the way it should and that it isn't disrupted by future additions to Rosetta.
// /// @author Sarah Biehn (biehn.4@osu.edu)

// // Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Package headers
#include <core/energy_methods/HRFDynamicsEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethod.hh>

//Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <basic/Tracer.hh>

//Utility headers
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>

static basic::Tracer TR("core.scoring.HRFDynamicsEnergyTests.cxxtest");

using namespace core;

class HRFDynamicsEnergyTests : public CxxTest::TestSuite {

public:

	HRFDynamicsEnergyTests() {};

	//Shared initialization
	void setUp() {
		core_init();
	}

	//Shared finalization
	void tearDown() {
	}

	//test score term
	void test_hrf_dynam() {

		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		TR << "Beginning HRFDynamicsEnergyTests::test_hrf_dynam() ..." << std::endl;

		core::scoring::methods::EnergyMethodOptions options;
		options.hrf_dynamics_input( "core/scoring/hrf_test.txt" ); //input HRF data to EMO options
		core::pose::Pose pose; //initialize pose
		core::import_pose::pose_from_file( pose, "core/scoring/hrf_dyn_test.pdb", core::import_pose::PDB_file ); //import PDB that corresponds to given HRF data
		ScoreFunction sfxn;
		sfxn.set_energy_method_options( options ); //give the score function the input file
		sfxn.set_weight( hrf_dynamics, 12.0 ); //rescore with hrf_dynamics
		TR << sfxn( pose ) << std::endl;
		//TR << "\t" << pose.sequence() << " \t " << 0 << " \t " << pose.energies() << std::endl;
		TR << "\t" << pose.energies().total_energies()[ hrf_dynamics ] << std::endl;
		TS_ASSERT_DELTA( sfxn( pose ), -198.358, 1e-1 ); //assess difference in scores
		TS_ASSERT_DELTA( pose.energies().total_energies()[ hrf_dynamics ], -16.5299, 1e-1 );
	}
};
