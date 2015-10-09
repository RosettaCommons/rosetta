// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/core/scoring/methods/AARepeatEnergy.cxxtest.hh
/// @brief  test suite for core::scoring::methods::AARepeatEnergy
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Test headers
#include <cxxtest/TestSuite.h>

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

#include <core/pose/annotated_sequence.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>

#include <core/pack/interaction_graph/ResidueArrayAnnealingEvaluator.hh>
#include <core/scoring/aa_repeat_energy/AARepeatEnergy.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

//Auto Headers
#include <utility/vector1.hh>

#include <map>

#include <boost/foreach.hpp>
#define foreach_         BOOST_FOREACH
#include <utility>


static basic::Tracer TR("core.scoring.methods.AARepeatEnergy.cxxtest");

// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::methods;
using namespace core::scoring::aa_repeat_energy;
using namespace core::scoring::annealing;

using namespace core::pack;
using namespace core::pack::task;
using namespace core::pack::rotamer_set;

class AARepeatEnergyTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void test_energy_annealing() {

		// Setup test pose
		Pose pose;
		make_pose_from_sequence( pose, "AAAAAAAA", "fa_standard");

		// Setup score function
		ScoreFunction scorefxn;
		scorefxn.set_weight( aa_repeat, 1 );
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
		for ( int r = 1; r < ev.get_num_nodes(); ++r ) {
			ev.set_state_for_node( r, r % 2 + 1);
		}
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 0, 1e-6);

		ev.set_state_for_node( 4, 2 );
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 1, 1e-6);

		PackerEnergy delta_energy;
		PackerEnergy pre_energy;

		ev.set_state_for_node( 4, 1 );
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 0, 1e-6);

		ev.consider_substitution( 4, 2, delta_energy, pre_energy );
		TS_ASSERT_DELTA( delta_energy, 1, 1e-6);
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 0, 1e-6);

		ev.consider_substitution( 6, 2, delta_energy, pre_energy );
		TS_ASSERT_DELTA( delta_energy, 1, 1e-6);
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 0, 1e-6);

		ev.commit_considered_substitution( );
		TS_ASSERT_DELTA( ev.get_energy_current_state_assignment(), 1, 1e-6);

		// Test via pack_rotamers run
		pack_rotamers(pose, scorefxn, task);
		PackerEnergy postpack_score = scorefxn(pose);

		TS_ASSERT_DELTA(postpack_score, 0, 1e-6);
	}

	void test_energy_cases() {
		typedef std::map< std::string, float> ScoreMap;
		ScoreMap expected_scores;

		expected_scores["A"] = 0;
		expected_scores["AA"] = 0;
		expected_scores["AAA"] = 1;
		expected_scores["AAAA"] = 10;
		expected_scores["AAAAA"] = 100;
		expected_scores["AAAAAA"] = 100;

		expected_scores["GAAAG"] = 1;
		expected_scores["GAAAAG"] = 10;
		expected_scores["GAAA"] = 1;
		expected_scores["AAAG"] = 1;

		expected_scores["AGAGAGAG"] = 0;
		expected_scores["AAAGAGAG"] = 1;
		expected_scores["AGAAAGAG"] = 1;
		expected_scores["AGAGAGGG"] = 1;

		expected_scores["GGGAAA"] = 2;
		expected_scores["GGGAAAA"] = 11;
		expected_scores["GGGGAAAA"] = 20;

		expected_scores["GGGGAAAATTTT"] = 30;
		expected_scores["GGGGAAAAGGGG"] = 30;
		expected_scores["GGGGAAAAGGGG"] = 30;

		ScoreFunction sfxn;
		sfxn.set_weight( aa_repeat, 1 );

		foreach_ ( ScoreMap::value_type &i, expected_scores ) {
			Pose test_pose;
			make_pose_from_sequence( test_pose, i.first, "fa_standard");

			TR << i.first << " " << i.second << std::endl;

			TS_ASSERT_DELTA( sfxn(test_pose), i.second, 1e-6);
		}
	}
	/// @brief Test the energy calculation for varying numbers of repeating residues.
	///
	void test_energy_eval() {
		if ( TR.visible() ) {
			TR << "Starting AARepeatEnergyTests::test_energy_eval()." << std::endl;
			TR << "Test created 20 July 2015 by Vikram K. Mulligan, Baker laboratory." << std::endl;
			TR << "This test checks that the aa_repeat_energy score term evaluates its energy correctly." << std::endl;
		}

		using namespace core::chemical;
		ResidueTypeSetCOP standard_residues( ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

		Pose trpcage( create_trpcage_ideal_pose() );
		ScoreFunction sfxn;
		sfxn.set_weight( aa_repeat, 0.5 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TEST\tEXPECTED\tACTUAL" << std::endl;
		if ( TR.visible() ) TR << "TrpCage:\t" << "0.5\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 0.5, 1e-6 );

		//Append three more alanines:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 20 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 20 );
		core::conformation::ResidueOP new_rsd1( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		core::conformation::ResidueOP new_rsd2( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		core::conformation::ResidueOP new_rsd3( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		trpcage.append_residue_by_bond(*new_rsd1, true, 0, 20, 0, false);
		trpcage.append_residue_by_bond(*new_rsd2, true, 0, 21, 0, false);
		trpcage.append_residue_by_bond(*new_rsd3, true, 0, 22, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 23 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 23 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+3:\t" << "1.0\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 1.0, 1e-6 );

		//Append one more alanine:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 23 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 23 );
		core::conformation::ResidueOP new_rsd4( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		trpcage.append_residue_by_bond(*new_rsd4, true, 0, 23, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 24 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 24 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+4:\t" << "5.5\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 5.5, 1e-6 );

		//Append one more alanine:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 24 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 24 );
		core::conformation::ResidueOP new_rsd5( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		trpcage.append_residue_by_bond(*new_rsd5, true, 0, 24, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 25 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 25 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+5:\t" << "50.5\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 50.5, 1e-6 );

		//Append one more alanine:
		core::pose::remove_variant_type_from_pose_residue( trpcage, CUTPOINT_UPPER, 25 );
		core::pose::remove_variant_type_from_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 25 );
		core::conformation::ResidueOP new_rsd6( core::conformation::ResidueFactory::create_residue( standard_residues->name_map("ALA") ) );
		trpcage.append_residue_by_bond(*new_rsd6, true, 0, 25, 0, false);
		core::pose::add_variant_type_to_pose_residue( trpcage, CUTPOINT_UPPER, 26 );
		core::pose::add_variant_type_to_pose_residue( trpcage, UPPER_TERMINUS_VARIANT, 26 );

		sfxn(trpcage);
		if ( TR.visible() ) TR << "TrpCage+6:\t" << "50.5\t" << trpcage.energies().total_energy() << std::endl;
		TS_ASSERT_DELTA( trpcage.energies().total_energy(), 50.5, 1e-6 );

		if ( TR.visible() ) {
			TR << "Test AARepeatEnergyTests::test_energy_eval() complete." << std::endl;
			TR.flush();
		}
		return;
	}

};


