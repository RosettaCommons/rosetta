// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/aa_composition_energy/AACompositionEnergy.cxxtest.hh
/// @brief  Test suite for core::scoring::aa_composition_energy::AACompositionEnergy, an energy term for controlling
/// sequence composition during design.
/// @details See also the core::conformation::symmetry::MirrorSymmetricConformation unit tests.  These have
/// another example of AAComposition being set up from code (with constraints attached to the pose).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <core/scoring/aa_composition_energy/AACompositionEnergySetup.hh>
#include <core/scoring/aa_composition_energy/AACompositionEnergy.hh>

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


static basic::Tracer TR("core.scoring.aa_composition_energy.AACompositionEnergy.cxxtest");

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

class AACompositionEnergyTests_linconsttail : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Test the tail functions with linear below, const above.
	///
	void test_tailfunctions_lin_const() {
		core_init_with_additional_options("-score:aa_composition_setup_file core/scoring/aa_composition_energy/tailfunction_linear.comp -out:levels core.scoring.aa_composition_energy.AACompositionEnergy:500");

		if ( TR.visible() ) {
			TR << "Starting test_tailfunctions_lin_const()." << std::endl;
		}

		// Set up score function
		ScoreFunction scorefxn;
		scorefxn.set_weight( aa_composition, 1 );

		{ //Scope 1: within test range
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 0, 1e-6);
		}

		{ //Scope 2: one below ideal, within test range
			Pose pose;
			make_pose_from_sequence( pose, "AAAAGGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 10, 1e-6);
		}

		{ //Scope 3: two below ideal, out of test range (linear region)
			Pose pose;
			make_pose_from_sequence( pose, "AAAGGGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 20, 1e-6);
		}

		{ //Scope 4: three below ideal, out of test range (linear region)
			Pose pose;
			make_pose_from_sequence( pose, "AAGGGGGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 30, 1e-6);
		}

		{ //Scope 5: one above ideal, within test range
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAAGGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 20, 1e-6);
		}

		{ //Scope 6: two above ideal, out of test range (constant region)
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAAAGGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 20, 1e-6);
		}

		{ //Scope 6: three above ideal, out of test range (constant region)
			Pose pose;
			make_pose_from_sequence( pose, "AAAAAAAAGG", "fa_standard");
			TS_ASSERT_DELTA(scorefxn(pose), 20, 1e-6);
		}

		if ( TR.visible() ) {
			TR << "Test test_tailfunctions_lin_const() complete." << std::endl;
			TR.flush();
		}

		return;
	}

};
