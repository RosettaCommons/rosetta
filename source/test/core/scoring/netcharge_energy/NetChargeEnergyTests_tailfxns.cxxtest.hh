// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/netcharge_energy/NetChargeEnergy.cxxtest.hh
/// @brief  Test suite for core::scoring::netcharge_energy::NetChargeEnergy, an energy term for controlling
/// net charge during design.
/// @details See also the core::conformation::symmetry::MirrorSymmetricConformation unit tests.  These have
/// another example of NetCharge being set up from code (with constraints attached to the pose).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <core/scoring/netcharge_energy/NetChargeEnergySetup.hh>
#include <core/scoring/netcharge_energy/NetChargeEnergy.hh>

// Unit headers

#include <platform/types.hh>

// Package Headers
#include <test/util/deriv_funcs.hh>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/pose/annotated_sequence.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pack/interaction_graph/ResidueArrayAnnealingEvaluator.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>


static basic::Tracer TR("core.scoring.netcharge_energy.NetChargeEnergyTests_tailfxns.cxxtest");

// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;
using namespace core::scoring::annealing;

using namespace core::pack;
using namespace core::pack::task;
using namespace core::pack::rotamer_set;

class NetChargeEnergyTests_tailfxns : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Test the energy calculation with linear tailfunctions.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_linear_tailfxns( ) {
		// Set up score function
		EnergyMethodOptions options;
		utility::vector1< std::string > files(1, "core/scoring/netcharge_energy/linear_tailfxns.charge");
		options.set_netcharge_setup_files(files);
		ScoreFunction scorefxn;
		scorefxn.set_energy_method_options(options);
		scorefxn.set_weight( netcharge, 1 );

		// Set up test pose
		Pose pose;
		make_pose_from_sequence( pose, "EEEEDDEEEE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 200.0, 1e-6 );
		make_pose_from_sequence( pose, "EEEEDDGEEE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 150.0, 1e-6 );
		make_pose_from_sequence( pose, "EEEEDDKEEE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 100.0, 1e-6 );
		make_pose_from_sequence( pose, "EEEEDDKGEE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 50.0, 1e-6 );
		make_pose_from_sequence( pose, "EEEEDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 60.0, 1e-6 );
		make_pose_from_sequence( pose, "EEGEDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 40.0, 1e-6 );
		make_pose_from_sequence( pose, "EEREDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 30.0, 1e-6 );
		make_pose_from_sequence( pose, "EERGDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 0.0, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 20.0, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKRGE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 25.0, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKRGG", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 30.0, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKRKG", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 35.0, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKRKR", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 40.0, 1e-6 );
	}

	/// @brief Test the energy calculation with constant tailfunctions.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_constant_tailfxns( ) {
		// Set up score function
		EnergyMethodOptions options;
		utility::vector1< std::string > files(1, "core/scoring/netcharge_energy/constant_tailfxns.charge");
		options.set_netcharge_setup_files(files);
		ScoreFunction scorefxn;
		scorefxn.set_energy_method_options(options);
		scorefxn.set_weight( netcharge, 1 );

		// Set up test pose
		Pose pose;
		make_pose_from_sequence( pose, "EEEEDDEEEE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 100.0, 1e-6 );
		make_pose_from_sequence( pose, "EEEEDDGEEE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 100.0, 1e-6 );
		make_pose_from_sequence( pose, "EEEEDDKEEE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 100.0, 1e-6 );
		make_pose_from_sequence( pose, "EEEEDDKGEE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 50.0, 1e-6 );
		make_pose_from_sequence( pose, "EEEEDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 60.0, 1e-6 );
		make_pose_from_sequence( pose, "EEGEDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 40.0, 1e-6 );
		make_pose_from_sequence( pose, "EEREDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 30.0, 1e-6 );
		make_pose_from_sequence( pose, "EERGDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 0.0, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 20.0, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKRGE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 25.0, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKRGG", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 25.0, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKRKG", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 25.0, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKRKR", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 25.0, 1e-6 );
	}

	/// @brief Test the energy calculation with quadratic tailfunctions.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	void test_quadratic_tailfxns( ) {
		// Set up score function
		EnergyMethodOptions options;
		utility::vector1< std::string > files(1, "core/scoring/netcharge_energy/quadratic_tailfxns.charge");
		options.set_netcharge_setup_files(files);
		ScoreFunction scorefxn;
		scorefxn.set_energy_method_options(options);
		scorefxn.set_weight( netcharge, 1 );

		// Set up test pose
		Pose pose;
		make_pose_from_sequence( pose, "EEEEDDEEEE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 233.33333333333333, 1e-6 );
		make_pose_from_sequence( pose, "EEEEDDGEEE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 161.11111111111111, 1e-6 );
		make_pose_from_sequence( pose, "EEEEDDKEEE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 100.0, 1e-6 );
		make_pose_from_sequence( pose, "EEEEDDKGEE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 50.0, 1e-6 );
		make_pose_from_sequence( pose, "EEEEDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 60.0, 1e-6 );
		make_pose_from_sequence( pose, "EEGEDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 40.0, 1e-6 );
		make_pose_from_sequence( pose, "EEREDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 30.0, 1e-6 );
		make_pose_from_sequence( pose, "EERGDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 0.0, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKREE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 20.0, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKRGE", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 25.0, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKRGG", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 33.33333333333333, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKRKG", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 45.0, 1e-6 );
		make_pose_from_sequence( pose, "EERKDDKRKR", "fa_standard");
		TS_ASSERT_DELTA( scorefxn(pose), 60.0, 1e-6 );
	}

};
