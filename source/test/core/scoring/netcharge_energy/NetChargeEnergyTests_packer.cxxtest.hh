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


static basic::Tracer TR("core.scoring.netcharge_energy.NetChargeEnergyTests_packer.cxxtest");

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

class NetChargeEnergyTests_packer : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options("-score:netcharge_setup_file negative_one.charge" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Test the energy calculation with the packer.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).  Modified from Alex Ford's tests for the aa_composition energy.
	void test_energy_annealing( ) {
		// Setup score function
		ScoreFunction scorefxn;
		scorefxn.set_weight( netcharge, 1 );

		// Setup test pose
		Pose pose;
		make_pose_from_sequence( pose, "DGGGEKGG", "fa_standard");
		TS_ASSERT_EQUALS(pose.residue_type(1).net_formal_charge(), -1);
		TS_ASSERT_EQUALS(pose.residue_type(2).net_formal_charge(),  0);
		TS_ASSERT_EQUALS(pose.residue_type(5).net_formal_charge(), -1);
		TS_ASSERT_EQUALS(pose.residue_type(6).net_formal_charge(),  1);
		TS_ASSERT_DELTA(scorefxn(pose), 0, 1e-6);

		make_pose_from_sequence( pose, "GGGGEKGG", "fa_standard");
		PackerEnergy prepack_score = scorefxn(pose);
		TS_ASSERT_DELTA(prepack_score, 10, 1e-6);

		make_pose_from_sequence( pose, "GGDDEKGG", "fa_standard");
		prepack_score = scorefxn(pose);
		TS_ASSERT_DELTA(prepack_score, 10, 1e-6);

		make_pose_from_sequence( pose, "GGDDEKEG", "fa_standard");
		prepack_score = scorefxn(pose);
		TS_ASSERT_DELTA(prepack_score, 40, 1e-6);

		make_pose_from_sequence( pose, "GGGGGKGG", "fa_standard");
		prepack_score = scorefxn(pose);
		TS_ASSERT_DELTA(prepack_score, 40, 1e-6);

		// Setup packer task and packer objects
		PackerTaskOP task( TaskFactory::create_packer_task( pose ));

		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		keep_aas[ core::chemical::aa_lys ] = true;
		keep_aas[ core::chemical::aa_asp ] = true;
		keep_aas[ core::chemical::aa_gly ] = true;
		keep_aas[ core::chemical::aa_glu ] = true;

		for ( core::Size i = 1; i <= pose.size(); i++ ) {
			task->nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas );
		}

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scorefxn, task );
		rotsets->build_rotamers( pose, scorefxn, packer_neighbor_graph );

		core::pack::interaction_graph::ResidueArrayAnnealingEvaluator ev;
		ev.initialize( scorefxn, pose, *rotsets, packer_neighbor_graph);

		// Base assignment should be base pose score...
		TS_ASSERT_EQUALS( ev.any_vertex_state_unassigned(), true );
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 40, 1e-6);

		for ( core::Size i(1); i<=5; ++i ) {
			// Test via pack_rotamers run five times.
			pack_rotamers(pose, scorefxn, task);
			PackerEnergy postpack_score = scorefxn(pose);
			TS_ASSERT_DELTA(postpack_score, 0, 1e-6);

			std::map< std::string, int > aa_count;
			aa_count["ASP"] = 0;
			aa_count["GLU"] = 0;
			aa_count["LYS"] = 0;
			aa_count["GLY"] = 0;

			for ( core::Size r = 1; r<= pose.size(); ++r ) {
				aa_count[ pose.residue(r).name3() ] += 1;
			}

			TS_ASSERT_EQUALS( aa_count["ASP"] + aa_count["GLU"] - aa_count["LYS"], 1);
		}
	}

};
