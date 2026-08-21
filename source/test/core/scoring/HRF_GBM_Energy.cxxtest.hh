// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
//
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// /// @file test/core/scoring/HRF_GBM_Energy.cxxtest.hh
// /// @brief Unit tests for HRF_GBM score term.
// /// @details This is a unit test for the HRF_GBM score function to ensure it's working as intended and that it isn't disrupted by future updates to Rosetta.
// /// @author Elijah Day (ehday@ucla.edu)

// // Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Package headers
#include <core/scoring/methods/EnergyMethodOptions.hh>

//Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <basic/Tracer.hh>

//Utility headers
#include <core/import_pose/import_pose.hh>

static basic::Tracer TR("core.scoring.HRF_GBM_EnergyTests.cxxtest");

using namespace core;

class HRF_GBM_EnergyTests : public CxxTest::TestSuite {

public:

	HRF_GBM_EnergyTests() {};

	//Shared initialization
	void setUp() {
		core_init();
	}

	//Shared finalization
	void tearDown() {
	}

	//test score term
	void test_hrf_gbm() {

		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		TR << "Beginning HRF_GBM_EnergyTests::test_hrf_gbm() ..." << std::endl;

		core::scoring::methods::EnergyMethodOptions options;
		options.hrf_gbm_input( "core/scoring/hrf_gbm_test.txt" ); //input sample HRF_GBM data to EMO options
		core::pose::Pose pose; //initialize pose
		core::import_pose::pose_from_file( pose, "core/scoring/hrf_gbm_test.pdb", core::import_pose::PDB_file ); //import PDB that corresponds to sample data
		ScoreFunction sfxn;
		sfxn.set_energy_method_options( options ); //pass the input file to the score function
		sfxn.set_weight( hrf_gbm, 9.0 ); //rescore with hrf_gbm
		TR << sfxn( pose ) << std::endl;
		TR << "\t" << pose.energies().total_energies()[ hrf_gbm ] << std::endl;
		TS_ASSERT_DELTA( sfxn( pose ), -74.818, 1e-1 ); //assess difference in scores
		TS_ASSERT_DELTA( pose.energies().total_energies()[ hrf_gbm ], -8.313, 1e-1 );
	}
};
