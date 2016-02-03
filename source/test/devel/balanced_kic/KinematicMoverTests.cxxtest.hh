// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/devel/balanced_kic/KinematicMover.cxxtest.hh
/// @brief  Weakly test the KinematicMover interface.
/// @author Kale Kundert
///
/// This test does not actually compare a move to a set of known results.  This 
/// kind of in-depth testing is already performed for each of the algorithm 
/// helper functions, and since the KinematicMover is just a composition of 
/// these algorithms, performing these tests again would be redundant.  
/// Instead, this test simply makes sure that all aspects of the mover 
/// interface can be invoked without crashing the program.

#ifndef INCLUDED_devel_balanced_kic_kinematic_mover_CXXTEST_HH
#define INCLUDED_devel_balanced_kic_kinematic_mover_CXXTEST_HH

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// The second include is neccesary to avoid an error about the incomplete 
// forward declaration of KinematicPerturber.  I'd like the KinematicMover 
// header to be self-sufficient, but I'm not sure how to accomplish that.
//
// I must actually be using the definition of a KinematicPerturber somewhere in 

#include <devel/balanced_kic/KinematicMover.hh>
#include <devel/balanced_kic/KinematicPerturber.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

using namespace core;
using namespace devel::balanced_kic;

class KinematicMoverTests : public CxxTest::TestSuite {

public:

	void setUp() { core_init(); }
	void tearDown() {}

	void test_mover_invocation() {
		pose::Pose pose;
		KinematicMover mover;

		import_pose::pose_from_file(pose, "devel/balanced_kic/loop.pdb", core::import_pose::PDB_file);

		mover.set_pivots(2, 3, 4);

		for (int i = 0; i < 1000; i++) {
			mover.apply(pose);
		}
	}

};

#endif


