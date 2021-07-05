// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
//
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// /// @file test/core/scoring/DEPC_MS_Energy.cxxtest.hh
// /// @brief Unit tests for depc_ms score term.
// /// @details This is a unit test for DEPC_MS_Energy (depc_ms score term) to ensure it's working the way it should and that it isn't disrupted by future additions to Rosetta.
// /// @author Sarah Biehn (biehn.4@osu.edu)

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

static basic::Tracer TR("core.scoring.DEPC_MS_EnergyTests.cxxtest");

using namespace core;

class DEPC_MS_EnergyTests : public CxxTest::TestSuite {

public:

	DEPC_MS_EnergyTests() {};

	//Shared initialization
	void setUp() {
		core_init();
	}

	//Shared finalization
	void tearDown() {
	}

	//test score term
	void test_depc_ms() {

		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		TR << "Beginning DEPC_MS_EnergyTests::test_depc_ms() ..." << std::endl;

		core::scoring::methods::EnergyMethodOptions options;
		options.depc_ms_input( "core/scoring/depc_test.txt" ); //input DEPC data to EMO options
		core::pose::Pose pose; //initialize pose
		core::import_pose::pose_from_file( pose, "core/scoring/hrf_dyn_test.pdb", core::import_pose::PDB_file ); //import PDB that corresponds to DEPC data; reuse the same PDB as HRF dynamics test because we have DEPC data for it too
		ScoreFunction sfxn;
		sfxn.set_energy_method_options( options ); //give the score function the input file
		sfxn.set_weight( depc_ms, 9.0 ); //rescore with depc_ms
		Real score( sfxn( pose ) );
		TR << score << std::endl;
		TR << "\t" << pose.energies().total_energies()[ depc_ms ] << std::endl;
		TS_ASSERT_DELTA( score, -116.717, 1e-1 ); //assess difference in scores
		TS_ASSERT_DELTA( pose.energies().total_energies()[ depc_ms ], -12.969, 1e-1 );
	}
};
