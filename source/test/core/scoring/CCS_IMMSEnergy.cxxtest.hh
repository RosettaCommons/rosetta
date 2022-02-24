// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
//
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// /// @file test/core/scoring/CCS_IMMSEnergy.cxxtest.hh
// /// @brief Unit tests for CCS_IMMSEnergy score term.
// /// @details This is a unit test for score function CCS_IMMSEnergy to ensure it's working the way it should and that it isn't disrupted by future additions to Rosetta.
// /// @author Bargeen Turzo (turzo.1@osu.edu)

// // Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Package headers
#include <core/energy_methods/CCS_IMMSEnergy.hh>
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

static basic::Tracer TR("core.scoring.CCS_IMMSEnergyTests.cxxtest");

using namespace core;

class CCS_IMMSEnergyTests : public CxxTest::TestSuite {

public:

	CCS_IMMSEnergyTests() {};

	//Shared initialization
	void setUp() {
		core_init_with_additional_options( "-ccs_exp 588 -ccs_nrots 300 -ccs_prad 1.0" );
	}

	//Shared finalization
	void tearDown() {
	}

	//test score term
	void test_ccs_imms() {

		using namespace core;
		using namespace core::pose;
		using namespace core::scoring;

		TR << "Beginning CCS_IMMSEnergyTests::test_ccs_imms() ..." << std::endl;
		core::pose::Pose pose; //initialize pose
		core::import_pose::pose_from_file( pose, "core/scoring/ccs_imms_test.pdb", core::import_pose::PDB_file ); //import PDB that corresponds to given CCS data
		ScoreFunction sfxn;
		sfxn.set_weight( ccs_imms, 1.0 ); //rescore with ccs_imms
		TR << sfxn( pose ) << std::endl;
		//TR << "\t" << pose.energies().total_energies()[ ccs_imms ] << std::endl;
		TS_ASSERT_DELTA( sfxn( pose ), 100, 1e-1 ); //assess difference in scores
	}
};
