// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/mhc_epitope_energy/MHCEpitopeEnergy.cxxtest.hh
/// @brief  Test suite for core::scoring::mhc_epitope_energy::MHCEptiopeEnergy, encapsulating prediction of MHC-peptide binding, to enable deimmunization by mutagenic epitope deletion
/// The code is largely based on (via copying and modifying) NetChargeEnergy (helpers) and HBNetEnergy (only updating around substitution position).
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

// Test headers
#include <cxxtest/TestSuite.h>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergy.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergyCreator.hh>

#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/protein_interface_design/movers/AddChainBreak.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/optimization/MinimizerMap.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/interaction_graph/ResidueArrayAnnealingEvaluator.hh>

// Package Headers
#include <test/core/init_util.hh>

#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>


static basic::Tracer TR("core.scoring.mhc_epitope_energy.MHCEpitopeEnergy.cxxtest");

// --------------- Test Class --------------- //

// using declarations
using namespace core;
using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::mhc_epitope_energy;

class MHCEpitopeEnergyTests_Scoring : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //
	// Shared variables
	Pose pose;
	core::Size nres;
	utility::vector1< core::conformation::ResidueCOP > reslist;
	core::Real seq_pp_score;

	// Shared initialization goes here.
	void setUp() {
		core_init();

		// Setup test pose
		std::string sequence = "YFCTRAFRILAWIGIQNPTS";
		make_pose_from_sequence( pose, sequence, "fa_standard");
		seq_pp_score = 11.0; // This is the propred score for the above sequence.

		// Setup a vector of of pointers from the pose.
		nres = pose.size();
		reslist.clear();
		reslist.reserve(nres);
		for ( core::Size ir=1; ir<=nres; ++ir ) {
			reslist.push_back( pose.residue(ir).get_self_ptr() );
		}
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Calculate the mhc_energy using a scorefunction and directly using calculate_energy.
	/// @brief Also test the disabling and re-enabling of the scoreterm with the minimizer.
	/// @author Brahm Yachnin
	void test_mhc_energy_matrix_scoring() {
		//// Calculate the energy by setting up an mhc_epitope scorefunction

		//Setup the config file
		utility::vector1< std::string > files(1, "propred8_5.mhc");
		methods::EnergyMethodOptions options;
		options.set_mhc_epitope_setup_files(files);

		//Associate the config with a scorefunction
		ScoreFunction scorefxn;
		scorefxn.set_weight( mhc_epitope, 1 );
		scorefxn.set_energy_method_options(options);

		// Calculate the score from the scorefxn.
		core::Real sfxn_energy;
		sfxn_energy = scorefxn(pose);
		TS_ASSERT_EQUALS( sfxn_energy, seq_pp_score );

		//// Now calculate the energy directly using MHCEpitopeEnergy's calculate_energy method

		// Set up a MHCEpitopeEnergy object
		MHCEpitopeEnergyOP mhc_energy( utility::pointer::make_shared<MHCEpitopeEnergy>( options ) );

		// Set up mhc_energy for packing
		core::pack::rotamer_set::RotamerSets rot_set; // This is needed for the following function, even though we don't use it.
		mhc_energy->set_up_residuearrayannealableenergy_for_packing( pose, rot_set, scorefxn );

		utility::vector1< core::Size > const dummy_rotamers( reslist.size(), 0 );

		// Calculate the score using the calculate_energy method.
		core::Real calc_energy;
		calc_energy = mhc_energy->calculate_energy( reslist, dummy_rotamers, 0 );
		TS_ASSERT_EQUALS( calc_energy, sfxn_energy );

		//// Disable in setup_for_minimizing and make sure the score is now 0.
		// Make a minmap to pass to setup_for_minimizing
		optimization::MinimizerMap minmap;
		// Disable the scoreterm for minimization
		mhc_energy->setup_for_minimizing( pose, scorefxn, minmap );
		// Calculate the score and make sure it's 0
		TS_ASSERT_EQUALS( mhc_energy->calculate_energy( reslist, dummy_rotamers, 0 ), 0 );

		// Re-enable post-minimization
		mhc_energy->finalize_after_minimizing( pose );
		TS_ASSERT_EQUALS( mhc_energy->calculate_energy( reslist, dummy_rotamers, 0 ), calc_energy );

		TR << "End of test_mhc_energy_matrix_scoring." << std::endl;
	}

	/// @brief Make a mutation, and check that the mhc_epitope cache is appropriately updated.
	/// @author Brahm Yachnin
	void test_mhc_energy_caching() {
		//// Make a substitution and update the score (using the caching system).

		//Setup the config file
		utility::vector1< std::string > files(1, "propred8_5.mhc");
		methods::EnergyMethodOptions options;
		options.set_mhc_epitope_setup_files(files);

		//Associate the config with a scorefunction
		ScoreFunction scorefxn;
		scorefxn.set_weight( mhc_epitope, 1 );
		scorefxn.set_energy_method_options(options);

		// Set up a MHCEpitopeEnergy object
		MHCEpitopeEnergyOP mhc_energy( utility::pointer::make_shared<MHCEpitopeEnergy>( options ) );

		utility::vector1< core::Size > const dummy_rotamers( reslist.size(), 0 );

		// Set up mhc_energy for packing and calculate the energy before mutation
		core::pack::rotamer_set::RotamerSets rot_set; // This is needed for the following function, even though we don't use it.
		mhc_energy->set_up_residuearrayannealableenergy_for_packing( pose, rot_set, scorefxn );
		mhc_energy->calculate_energy( reslist, dummy_rotamers, 0 );

		// The position and new ID of the substituted residue.
		core::Size subst_position = 9;
		std::string subst_id = "GLU";
		// Make and apply the MutateResidue mover
		protocols::simple_moves::MutateResidue mutate_res( subst_position, subst_id );
		mutate_res.apply(pose);

		// Reset the reslist vector
		reslist.clear();
		reslist.reserve(nres);
		for ( core::Size ir=1; ir<=nres; ++ir ) {
			reslist.push_back( pose.residue(ir).get_self_ptr() );
		}

		core::Real cache_energy, full_energy;
		// Calculate the energy by updating the cache at the substitution position only
		cache_energy = mhc_energy->calculate_energy( reslist, dummy_rotamers, subst_position );
		TS_ASSERT_EQUALS( cache_energy, 7.0 );
		// Calculate the energy of the full pose, and verify that it's the same
		full_energy = mhc_energy->calculate_energy( reslist, dummy_rotamers, 0 );
		TS_ASSERT_EQUALS( full_energy, cache_energy );

		TR << "End of test_mhc_energy_caching." << std::endl;
	}

	/// @brief Check the function of the external database, including mutations.
	/// @author Brahm Yachnin
	void test_mhc_energy_external() {
#ifndef MULTI_THREADED
		//Setup the config file
		utility::vector1< std::string > files(1, "core/scoring/mhc_epitope_energy/external_db.mhc");
		methods::EnergyMethodOptions options;
		options.set_mhc_epitope_setup_files(files);

		//Associate the config with a scorefunction
		ScoreFunction scorefxn;
		scorefxn.set_weight( mhc_epitope, 1 );
		scorefxn.set_energy_method_options(options);

		// Set up a MHCEpitopeEnergy object
		MHCEpitopeEnergyOP mhc_energy( utility::pointer::make_shared<MHCEpitopeEnergy>( options ) );

		// Score the pose
		TS_ASSERT_EQUALS(scorefxn(pose), 15);

		// Mutate using the MutateResidue mover
		// The position and new ID of the substituted residue.
		core::Size subst_position = 9;
		std::string subst_id = "ASP";
		// Make and apply the MutateResidue mover
		protocols::simple_moves::MutateResidue mutate_res( subst_position, subst_id );
		mutate_res.apply(pose);

		// Score the pose
		TS_ASSERT_EQUALS(scorefxn(pose), 8);

		// Mutate using the MutateResidue mover to an unknown residue
		// The position and new ID of the substituted residue.
		// subst_position = 9;
		subst_id = "ALA";
		// Make and apply the MutateResidue mover
		mutate_res.set_res_name(subst_id);
		mutate_res.apply(pose);

		// Score the pose
		TS_ASSERT_EQUALS(scorefxn(pose), 600);

#endif
		TR << "End of test_mhc_energy_external." << std::endl;
	}

	/// @brief Compare the scores for a pose before and after a chainbreak is introduced.
	/// @author Brahm Yachnin
	void test_mhc_energy_chainbreak() {
		//Setup the config file
		utility::vector1< std::string > files(1, "propred8_5.mhc");
		methods::EnergyMethodOptions options;
		options.set_mhc_epitope_setup_files(files);

		//Associate the config with a scorefunction
		ScoreFunction scorefxn;
		scorefxn.set_weight( mhc_epitope, 1 );
		scorefxn.set_energy_method_options(options);

		//Calculate the score from the scorefxn.
		core::Real fullpose_energy;
		fullpose_energy = scorefxn(pose);
		TS_ASSERT_EQUALS( fullpose_energy, seq_pp_score );

		//Store the sequence
		std::string fullpose_sequence = pose.sequence();

		//Introduce a chainbreak in pose at position break_pos;
		//Setup a mover to break the pose at position 10, and change the conformation
		std::string break_pos = "10";
		protocols::protein_interface_design::movers::AddChainBreak add_chainbreak;
		add_chainbreak.resnum(break_pos);
		add_chainbreak.change_conformation(true);
		//Apply the mover to pose
		add_chainbreak.apply(pose);

		//Assert that the full pose sequence is the same as the two by chain sequences combined
		TS_ASSERT_EQUALS( pose.chain_sequence(1)+pose.chain_sequence(2), fullpose_sequence );

		//Create a vector of poses and put each chain in its own pose.
		utility::vector1 < PoseOP > poses_by_chain;
		poses_by_chain = pose.split_by_chain();

		//Calculate the energy for the broken pose.
		core::Real broken_pose_energy = scorefxn(pose);
		//The chainbreak should disrupt epitopes that spanned the new chainbreak, so should be less than the full pose.
		TS_ASSERT_LESS_THAN(broken_pose_energy, fullpose_energy);

		//Calculate the energy for the individual chain poses, and store in energy_by_chain.
		utility::vector1 < core::Real > energy_by_chain;
		core::Real sum_energy_by_chain = 0;
		for ( core::Size i = 1; i <= poses_by_chain.size(); ++i ) {
			energy_by_chain.push_back( scorefxn(*poses_by_chain[i]) );
			sum_energy_by_chain += energy_by_chain[i];
		}

		//The broken pose should have scores equal to the individual chain poses.
		TS_ASSERT_EQUALS(broken_pose_energy, sum_energy_by_chain);
		TS_ASSERT_EQUALS(energy_by_chain[1], 2.0);
		TS_ASSERT_EQUALS(energy_by_chain[2], 3.0);

		TR << "End of test_mhc_energy_chainbreak." << std::endl;
	}

	/// @brief Compare the scores for a pose before and after a unnatural AA is introduced.
	/// @author Brahm Yachnin
	void test_mhc_energy_uaa() {
		//Setup the config file
		utility::vector1< std::string > files(1, "propred8_5.mhc");
		methods::EnergyMethodOptions options;
		options.set_mhc_epitope_setup_files(files);

		//Associate the config with a scorefunction
		ScoreFunction scorefxn;
		scorefxn.set_weight( mhc_epitope, 1 );
		scorefxn.set_energy_method_options(options);

		//Calculate the score from the scorefxn.
		core::Real fullpose_energy;
		fullpose_energy = scorefxn(pose);
		TS_ASSERT_EQUALS( fullpose_energy, seq_pp_score );

		// The position and new ID of the substituted residue.
		core::Size subst_position = 10;
		// B3L is the beta-3 version of Leu, the AA previously at position 10.
		std::string subst_id = "C95";
		// Make and apply the MutateResidue mover
		protocols::simple_moves::MutateResidue mutate_res( subst_position, subst_id );
		mutate_res.apply(pose);

		//Calculate the energy for the UAA pose.
		core::Real uaa_pose_energy = scorefxn(pose);
		//The UAA should disrupt epitopes that spanned it position, so score should be less than the full pose.
		TS_ASSERT_LESS_THAN(uaa_pose_energy, fullpose_energy);
		TS_ASSERT_EQUALS(uaa_pose_energy, 5.0);

		// Create new poses that contain only the residues before and after the UAA, respectively.
		Pose N_of_uaa = pose;
		Pose C_of_uaa = pose;
		N_of_uaa.delete_residue_range_slow( subst_position, N_of_uaa.size() );
		C_of_uaa.delete_residue_range_slow( 1, subst_position );

		// The uaa should break any epitopes, so the uaa_pose_energy should equal the sum of N+C_of_uaa pose energies.
		core::Real pose_fragment_energy_sum = scorefxn(N_of_uaa) + scorefxn(C_of_uaa);
		TS_ASSERT_EQUALS(uaa_pose_energy, pose_fragment_energy_sum);

		TR << "End of test_mhc_energy_uaa." << std::endl;
	}

	/// @brief Make the pose symmetric, and check that the score is properly updated.
	/// @author Brahm Yachnin
	void test_mhc_energy_sym() {
		//Setup the config file
		utility::vector1< std::string > files(1, "propred8_5.mhc");
		methods::EnergyMethodOptions options;
		options.set_mhc_epitope_setup_files(files);

		//Associate the config with the scorefunction
		ScoreFunction scorefxn;
		scorefxn.set_weight( mhc_epitope, 1 );
		scorefxn.set_energy_method_options(options);

		//Set up a MHCEpitopeEnergy object
		MHCEpitopeEnergyOP mhc_energy( utility::pointer::make_shared<MHCEpitopeEnergy>( options ) );
		//Set up mhc_energy for packing and calculate the energy before symmetrizing.
		core::pack::rotamer_set::RotamerSets rot_set; // This is needed for the following function, even though we don't use it.
		mhc_energy->set_up_residuearrayannealableenergy_for_packing( pose, rot_set, scorefxn );
		utility::vector1< core::Size > const dummy_rotamers( reslist.size(), 0 );
		core::Real presym_score = mhc_energy->calculate_energy( reslist, dummy_rotamers, 0 );

		//Make a symmetric pose, using a de novo c3 symmetry file.
		core::Size nsub = 3; //Three subunits in the c3.sym file.
		Pose sym_pose = pose;
		core::pose::symmetry::make_symmetric_pose( sym_pose, "core/scoring/mhc_epitope_energy/c3.sym" );
		//Score the pose
		core::Real sym_score = scorefxn(sym_pose);

		//The score should be the pre-sym score * number of symmetric subunits.
		TS_ASSERT_EQUALS( sym_score, presym_score * nsub );

		TR << "End of test_mhc_energy_sym." << std::endl;
	}

	/// @brief Test the mhc_energy with the packer.
	/// @details Largely copied from AACompositionEnergy_packer.cxxtest.hh
	/// @author Brahm Yachnin
	void test_mhc_energy_packer() {
		//Namespace
		using namespace core::pack;
		using namespace core::pack::task;
		using namespace core::pack::rotamer_set;

		//Setup the config file
		utility::vector1< std::string > files(1, "propred8_5.mhc");
		methods::EnergyMethodOptions options;
		options.set_mhc_epitope_setup_files(files);

		//Associate the config with the scorefunction
		ScoreFunction scorefxn;
		scorefxn.set_weight( mhc_epitope, 1 );
		scorefxn.set_energy_method_options(options);

		//Calculate the energy of the starting pose
		PackerEnergy prepack_score = scorefxn(pose);
		TS_ASSERT_EQUALS( prepack_score, seq_pp_score ); // Same as in previous tests.

		// Setup packer task
		PackerTaskOP task( TaskFactory::create_packer_task( pose ));

		// Allow the packer to sample the current residue, or ASP
		// ASP is typically not present in many epitopes, so should allow score to be optimized
		utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, false );
		keep_aas[ core::chemical::aa_asp ] = true;

		// Update the packer task to use only the native AA or ASP.
		// Also keep track of the number of initial ASPs in the pose.
		core::Size num_init_asp = 0;
		for ( core::Size i = 1; i <= pose.size(); i++ ) {
			task->nonconst_residue_task( i ).restrict_nonnative_canonical_aas( keep_aas );
			if ( pose.aa(i) == core::chemical::aa_asp ) {
				++num_init_asp;
			}
		}

		// Setup rotamer sets
		RotamerSetsOP rotsets( utility::pointer::make_shared<RotamerSets>() );
		rotsets->set_task( task );
		utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scorefxn, task );
		/* Ooops, get_n_rotamers_for_residue_type is from class RotamerSet, not RotamerSets.
		Would be nice if there was a more elegant way of doing this, rather than hardcoding the number of rotamers below.
		// Count the number of rotamers allowed for each position and store in nrots.
		// Get the nrots for the native, and also for ASP at each position
		core::Size nrots = 0;
		for ( core::Size i = 1; i <= pose.size(); i++ ) {
		nrots += rotsets->get_n_rotamers_for_residue_type( pose.aa(i) );
		if ( pose.aa(i) != core::chemical::aa_asp ) {
		nrots += rotsets->get_n_rotamers_for_residue_type( core::chemical::aa_asp );
		}
		}
		TR << "Calculated nrots: " << nrots << std::endl;
		*/
		rotsets->build_rotamers( pose, scorefxn, packer_neighbor_graph );

		// Figure out how many rotamers
		core::Size nrots = rotsets->nrotamers();
		// TS_ASSERT_EQUALS( nrots, 655 ); // This test is unstable depending on the compiler, compile mode, etc., so disabling.  (The problem exists without invoking any mhc_energy-related code.)

		// Initialize the interaction graph
		core::pack::interaction_graph::ResidueArrayAnnealingEvaluator ev;
		ev.initialize( scorefxn, pose, *rotsets, packer_neighbor_graph);

		//Data about the ev
		// The number of nodes (equals the number of residues in the pose)
		TS_ASSERT_EQUALS( ev.get_num_nodes(), pose.size() );
		// The number of states (equals the number of rotamers across all residues)
		TS_ASSERT_EQUALS( ev.get_num_total_states(), nrots );
		// The number of states (rotamers) for a particular node/residue, say 1
		TS_ASSERT_EQUALS( ev.get_num_states_for_node(1), 28 );
		// Check if any vertex state is unassigned.  We haven't assigned any of them, so should be true
		TS_ASSERT_EQUALS( ev.any_vertex_state_unassigned(), 1 );
		// Get the energy from the AnnealingEvaluator.  Since we haven't done anything, it should still equal prepack_score
		TS_ASSERT_EQUALS( ev.get_energy_current_state_assignment(), prepack_score );

		// Test state assignment and consideration
		// Loop over all positions, and mutate to D.
		for ( int r = 1; r <= ev.get_num_nodes(); ++r ) {
			//The list of states is based on alphabetically order of residues (right?)
			//If we want poly-D, all residues should be changed to state=1 unless it's ALA or CYS
			//For ALA and CYS, state should be the max.
			if ( pose.aa(r) == core::chemical::aa_ala || pose.aa(r) == core::chemical::aa_cys ) {
				ev.set_state_for_node( r, ev.get_num_states_for_node(r));
			} else {
				ev.set_state_for_node( r, 1 );
			}
		}

		// We've changed to poly-D.  This should have an mhc score of 0.
		TS_ASSERT_EQUALS( ev.get_energy_current_state_assignment(), 0 );

		// Loop over all positions, mutating back to the original sequence.
		for ( int r = 1; r <= ev.get_num_nodes(); ++r ) {
			//The list of states is based on alphabetically order of residues (right?)
			//Basically, just do the reverse of the above for loop.
			if ( pose.aa(r) == core::chemical::aa_ala || pose.aa(r) == core::chemical::aa_cys ) {
				ev.set_state_for_node( r, 1 );
			} else {
				ev.set_state_for_node( r, ev.get_num_states_for_node(r) );
			}
		}
		// We've changed back to the original sequence, so we should be back to prepack_score.
		TS_ASSERT_EQUALS( ev.get_energy_current_state_assignment(), prepack_score );

		// Values for the delta_energy and pre_energy, which are needed for consider_substitution.
		PackerEnergy delta_energy;
		PackerEnergy pre_energy;

		// Consider mutation residue 10 (Leu) to Asp.
		ev.consider_substitution(10, 1, delta_energy, pre_energy);

		// The mutation is not yet applied, so current energy should be unchanged
		TS_ASSERT_EQUALS( ev.get_energy_current_state_assignment(), prepack_score );

		// Commit the change.  The energy should now be prepack_score + delta_energy
		ev.commit_considered_substitution();
		TS_ASSERT_EQUALS( ev.get_energy_current_state_assignment(), prepack_score + delta_energy );

		// Test via pack_rotamers run
		pack_rotamers(pose, scorefxn, task);
		PackerEnergy postpack_score = scorefxn(pose);
		// The packer should be able to reduce the score to 0.
		TS_ASSERT_EQUALS( postpack_score, 0 );

		TR << "End of test_mhc_energy_packer." << std::endl;
	}

};
