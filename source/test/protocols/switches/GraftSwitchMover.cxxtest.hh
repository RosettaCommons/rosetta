// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/switches/GraftSwitchMover.cxxtest.hh
/// @brief unit tests for GraftSwitchMover
/// @author Bobby Langan (robert.langan@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <util/pose_funcs.hh>

// Project headers
#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>

//core
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

//class headers
#include <protocols/switches/GraftSwitchMover.hh>

static basic::Tracer TR("protocols.GraftSwitchMover.cxxtest.hh");

namespace {
class GraftSwitchMoverTests : public CxxTest::TestSuite {

public:

	// Shared set up goes here
	void setUp()
	{
		core_init();
	}

	void test_graft_switch_mover_setup()
	{
		using namespace protocols::switches;

		TR << "Start Test 1" << std::endl;
		TR << "Importing Pose" << std::endl;
		core::pose::Pose testPose;
		core::import_pose::pose_from_file( testPose, "protocols/switches/GraftSwitchMoverTest.pdb", core::import_pose::PDB_file);

		//Create the mover
		TR << "Creating mover" << std::endl;
		//HBNet hbnet_mover;
		GraftSwitchMover lockr_mover;
		TR << "Intializing mover" << std::endl;
		//std::set< Size > const start{46};
		//hbnet_mover.set_start_resnums( start );
		lockr_mover.add_sequence("L-MSCAQES");
		lockr_mover.set_start_graft(191);
		lockr_mover.set_end_graft(199);
		utility::vector1< core::Size > i_res;
		i_res.push_back(5);
		i_res.push_back(6);
		lockr_mover.set_important_residues(i_res);
		lockr_mover.set_pack_neighbors(false);
		lockr_mover.set_pack_min(false);
		lockr_mover.score_function( core::scoring::get_score_function() );

		TR << "Applying mover" << std::endl;
		//Apply the mover, will call all key functions
		lockr_mover.apply(testPose);

		//all of most important functions called by .apply() above, but should add more individual funciton tests here:

		//checks
		//check that we found a thread mover was applied correctly
		TS_ASSERT( testPose.sequence().find("MSCAQES") > 190 );

		TR << "End Test" << std::endl;
	}

};

}

