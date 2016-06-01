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

class AACompositionEnergyTests_packer : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Test the energy calculation with the packer.
	/// @author Alex Ford.
	void test_energy_annealing( ) {
		core_init_with_additional_options("-score:aa_composition_setup_file exactly_one_ala.comp");

		// Setup score function
		ScoreFunction scorefxn;
		scorefxn.set_weight( aa_composition, 1 );

		// Setup test pose
		Pose pose;
		make_pose_from_sequence( pose, "AGGGGGGG", "fa_standard");
		TS_ASSERT_DELTA(scorefxn(pose), 0, 1e-6);

		make_pose_from_sequence( pose, "AAAAAAAA", "fa_standard");

		PackerEnergy prepack_score = scorefxn(pose);
		TS_ASSERT_DELTA(prepack_score, 100, 1e-6);

		// Setup packer task and packer objects
		PackerTaskOP task( TaskFactory::create_packer_task( pose ));

		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		keep_aas[ core::chemical::aa_ala ] = true;
		keep_aas[ core::chemical::aa_gly ] = true;

		for ( core::Size i = 1; i <= pose.total_residue(); i++ ) {
			task->nonconst_residue_task( i ).restrict_absent_canonical_aas( keep_aas );
		}

		RotamerSetsOP rotsets( new RotamerSets() );
		rotsets->set_task( task );
		graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scorefxn, task );
		rotsets->build_rotamers( pose, scorefxn, packer_neighbor_graph );

		TS_ASSERT_EQUALS( rotsets->nrotamers(), pose.total_residue() * 2);

		core::pack::interaction_graph::ResidueArrayAnnealingEvaluator ev;
		ev.initialize( scorefxn, pose, *rotsets, packer_neighbor_graph);

		TS_ASSERT_EQUALS( ev.get_num_nodes(), 8);
		TS_ASSERT_EQUALS( ev.get_num_total_states(), 16);
		TS_ASSERT_EQUALS( ev.get_num_states_for_node(1), 2);

		// Base assignment should be base pose score...
		TS_ASSERT_EQUALS( ev.any_vertex_state_unassigned(), true );
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 100, 1e-6);

		// Test state assignment and consideration
		ev.set_state_for_node(1, 1);
		for ( int r = 2; r <= ev.get_num_nodes(); ++r ) {
			ev.set_state_for_node( r, 2);
		}

		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 0, 1e-6);

		ev.set_state_for_node(2, 1);
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 100, 1e-6);

		PackerEnergy delta_energy;
		PackerEnergy pre_energy;

		ev.consider_substitution(2, 2, delta_energy, pre_energy);

		TS_ASSERT_DELTA( delta_energy, -100, 1e-6);
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 100, 1e-6);

		ev.commit_considered_substitution();
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 0, 1e-6);

		// Test via pack_rotamers run
		pack_rotamers(pose, scorefxn, task);
		PackerEnergy postpack_score = scorefxn(pose);
		TS_ASSERT_DELTA(postpack_score, 0, 1e-6);

		std::map< std::string, int > aa_count;
		aa_count["ALA"] = 0;
		aa_count["GLY"] = 0;

		for ( core::Size r = 1; r<= pose.total_residue(); ++r ) {
			aa_count[ pose.residue(r).name3() ] += 1;
		}

		TS_ASSERT_EQUALS( aa_count["ALA"], 1);
		TS_ASSERT_EQUALS( aa_count["GLY"], pose.total_residue() - 1);
	}

};
