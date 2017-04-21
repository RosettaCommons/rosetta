// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/hbnet/HBNet.cxxtest.hh
/// @brief unit tests for HBNet
/// @author Scott Boyken (sboyken@gmail.com)

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

//class headers
#include <protocols/hbnet/HBNet.hh>
#include <protocols/hbnet/HBNetStapleInterface.hh>
//#include <protocols/hbnet/HBNet_util.hh>
//#include <protocols/simple_moves/symmetry/DetectSymmetryMover.hh>

static basic::Tracer TR("protocols.HBNet.cxxtest.hh");

namespace {
class HBNetMoverTests : public CxxTest::TestSuite {

public:

	// Shared set up goes here
	void setUp()
	{
		core_init();
	}

	void test_hbnet_mover_setup()
	{
		using namespace protocols::hbnet;

		TR << "Start Test 1" << std::endl;

		//Import the test scaffold
		TR << "Importing Pose" << std::endl;
		core::pose::Pose testPose;
		core::import_pose::pose_from_file( testPose, "protocols/hbnet/SB13.pdb" , core::import_pose::PDB_file);

		//protocols::simple_moves::symmetry::DetectSymmetry detect_symm;
		//detect_symm.apply(testPose);

		//Create the mover
		TR << "Creating mover" << std::endl;
		//HBNet hbnet_mover;
		HBNetStapleInterface hbnet_mover;
		TR << "Intializing mover" << std::endl;
		//std::set< Size > const start{46};
		//hbnet_mover.set_start_resnums( start );
		hbnet_mover.set_find_only_native(true);
		hbnet_mover.set_min_networks_size( 4 );
		hbnet_mover.set_max_unsat( 3 );

		TR << "Applying mover" << std::endl;
		//Apply the mover, will call all key functions
		hbnet_mover.apply(testPose);

		//all of most important functions called by .apply() above, but should add more individual funciton tests here:

		//checks
		//check that we found netowrks and mover was applied correctly
		TS_ASSERT( hbnet_mover.get_native_vec().size() > 0 );

		TR << "End Test" << std::endl;
	}
};
}

