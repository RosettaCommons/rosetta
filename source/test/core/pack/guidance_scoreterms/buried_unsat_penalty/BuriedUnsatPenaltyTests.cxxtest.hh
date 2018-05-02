// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pack/guidance_scoreterms/buried_unsat_penalty/BuriedUnsatPenaltyTests.cxxtest.hh
/// @brief  Unit tests for the BuriedUnsatPenalty energy method, which provides a packer-compatible penalty for buried unsatisfied Hbond donors or acceptors.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/BuriedUnsatPenalty.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraph.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphContainer.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/ResidueLevelTask.hh>
#include <core/pack/task/ResidueLevelTask_.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

// Protocols Headers, for convenience
#include <protocols/simple_moves/MutateResidue.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("BuriedUnsatPenaltyTests");


class BuriedUnsatPenaltyTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}

	/// @brief Simple test: just score a pose.
	/// @details This pose was created with Scott's HBNet, and should have fully satisfied hydrogen bonds.
	void test_simple_scoring(){
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/guidance_scoreterms/buried_unsat_penalty/graph/hbnet_testcase1.pdb"  );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::buried_unsatisfied_penalty, 1.0 );

		try {
			sfxn(pose);
		} catch( utility::excn::Exception & exception ) {
			TR << "Exception caught in unit test: " << exception.msg() << std::endl;
			TS_ASSERT(false);
		} catch ( ... ) {
			TR << "Unhandled exception caught in unit test!" << std::endl;
			TS_ASSERT(false);
		}

		core::Real const EXPECTED_ENERGY( std::pow(9*5, 2) ); //Nine unsats.
		TS_ASSERT_DELTA( pose.energies().total_energy(), EXPECTED_ENERGY, 1e-6 );

		protocols::simple_moves::MutateResidue mutres( 54, "ALA" ); //Mutate a hydrogen bond-forming residue.
		mutres.apply(pose);
		// This should result in 2 acceptors unsatisfied.  1 donor group loses a hydrogen bond, but that's an "acceptable" loss.
		sfxn(pose); //Rescore
		TS_ASSERT_DELTA( pose.energies().total_energy(), static_cast< core::Real >( std::pow( 11*5, 2 ) ), 1e-6 );
		//pose.dump_pdb( "BURIED_UNSAT_TEST_SIMPLE_SCORING.pdb" ); //DELETE ME

		//The following nonsense is just to dump out the groups:
		core::scoring::methods::EnergyMethodOptions options;
		core::pack::guidance_scoreterms::buried_unsat_penalty::BuriedUnsatPenalty scoreterm( options );
		core::pack::task::PackerTaskOP task( new core::pack::task::PackerTask_(pose) );
		for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
			(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(i))).restrict_to_repacking();
		}
		core::pack::rotamer_set::RotamerSetsOP rotsets( utility::pointer::make_shared< core::pack::rotamer_set::RotamerSets >() );
		rotsets->set_task( task );
		rotsets->initialize_pose_for_rotsets_creation(pose);
		utility::graph::GraphOP packer_neighbor_graph( core::pack::create_packer_graph( pose, sfxn, task ) );
		rotsets->build_rotamers(pose, sfxn, packer_neighbor_graph);
		rotsets->prepare_sets_for_packing(pose, sfxn);
		scoreterm.set_up_residuearrayannealableenergy_for_packing( pose, *rotsets, sfxn );
		scoreterm.provide_pymol_commands_to_show_groups( TR, pose );
		TR.flush();
	}

	/// @brief Test penalty graph setup for packing.
	void test_setup_for_packing() {
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/guidance_scoreterms/buried_unsat_penalty/graph/hbnet_testcase1.pdb"  );

		core::scoring::methods::EnergyMethodOptions options;
		core::pack::guidance_scoreterms::buried_unsat_penalty::BuriedUnsatPenalty scoreterm( options );
		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight(core::scoring::buried_unsatisfied_penalty, 1.0);
		//sfxn->set_weight( core::scoring::hbond_sc, 1.0 );

		(*sfxn)(pose);

		core::pack::task::PackerTaskOP task( new core::pack::task::PackerTask_(pose) );
		for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
			(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(i))).restrict_to_repacking();
		}

		core::pack::rotamer_set::RotamerSetsOP rotsets( new core::pack::rotamer_set::RotamerSets );
		rotsets->set_task( task );
		rotsets->initialize_pose_for_rotsets_creation(pose);
		utility::graph::GraphOP packer_neighbor_graph( core::pack::create_packer_graph( pose, *sfxn, task ) );
		rotsets->build_rotamers(pose, *sfxn, packer_neighbor_graph);
		rotsets->prepare_sets_for_packing(pose, *sfxn);

		scoreterm.set_up_residuearrayannealableenergy_for_packing( pose, *rotsets, *sfxn );
		core::pack::guidance_scoreterms::buried_unsat_penalty::graph::BuriedUnsatPenaltyGraph const & bunsat_graph( *(scoreterm.unsat_graph_) );
		TS_ASSERT_EQUALS( rotsets->total_residue(), pose.total_residue() ); //Should be true
		for ( core::Size i(1), imax(rotsets->total_residue()); i<=imax; ++i ) {
			core::pack::rotamer_set::RotamerSetCOP rotset( rotsets->rotamer_set_for_residue(i) );
			TR << "Residue " << i << " has " <<  rotset->num_rotamers() << " rotamers." << std::endl;
			for ( core::Size j(1), jmax(rotset->num_rotamers()); j<=jmax; ++j ) {
				core::conformation::ResidueCOP rotamer( rotset->rotamer(j) );
				TS_ASSERT( rotamer != nullptr );
				core::Size const nodeindex( bunsat_graph.get_node_index( rotamer ) );
				TS_ASSERT( nodeindex != 0 ); //Ensure that each rotamer has a corresponding node in the graph.
				TS_ASSERT( bunsat_graph.residue_memory_address_to_nodeindex_.at( rotamer ) == nodeindex );
				TS_ASSERT( bunsat_graph.nodeindex_to_residue_memory_address( nodeindex ) == rotamer ); //Ensure that we're storing the memory address/node index association
			}
			TS_ASSERT( bunsat_graph.residue_memory_address_to_nodeindex_.at( pose.conformation().residue_cop(i) ) != 0 );
		}

		scoreterm.provide_pymol_commands_to_show_groups( TR, pose );
		//pose.dump_pdb( "BURIED_UNSAT_TESTPOSE_SET_UP_FOR_PACKING.pdb" ); //DELETE ME.
		TR << std::endl;
		TR.flush();
	}

	/// @brief Test the function calls made during a packing trajectory.
	void test_packing_steps() {
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/guidance_scoreterms/buried_unsat_penalty/graph/hbnet_testcase2.pdb"  );

		core::scoring::methods::EnergyMethodOptions options;
		options.buried_unsatisfied_penalty_burial_threshold(1.0);
		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::buried_unsatisfied_penalty, 1.0 );
		//sfxn->set_weight( core::scoring::hbond_sc, 1.0 );
		sfxn->set_energy_method_options(options);

		(*sfxn)(pose);

		core::pack::guidance_scoreterms::buried_unsat_penalty::BuriedUnsatPenalty scoreterm(options), scoreterm2(options);

		core::pack::task::PackerTaskOP task( new core::pack::task::PackerTask_(pose));

		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(1))).restrict_to_repacking();
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(3))).restrict_to_repacking();
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(4))).restrict_to_repacking();
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(6))).restrict_to_repacking();
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(7))).restrict_to_repacking();
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(9))).restrict_to_repacking();

		utility::vector1< bool > allowed_aas( core::chemical::num_canonical_aas, false );
		allowed_aas[ core::chemical::aa_asn ] = true;
		allowed_aas[ core::chemical::aa_asp ] = true;
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(2))).restrict_absent_canonical_aas( allowed_aas );
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(5))).restrict_absent_canonical_aas( allowed_aas );
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(8))).restrict_absent_canonical_aas( allowed_aas );

		core::pack::rotamer_set::RotamerSetsOP rotsets( new core::pack::rotamer_set::RotamerSets );
		rotsets->set_task( task );
		rotsets->initialize_pose_for_rotsets_creation(pose);
		utility::graph::GraphOP packer_neighbor_graph( core::pack::create_packer_graph( pose, *sfxn, task ) );
		rotsets->build_rotamers(pose, *sfxn, packer_neighbor_graph);
		rotsets->prepare_sets_for_packing(pose, *sfxn);

		scoreterm.set_up_residuearrayannealableenergy_for_packing( pose, *rotsets, *sfxn );
		scoreterm.provide_pymol_commands_to_show_groups( TR, pose );
		TR.flush();
		//pose.dump_pdb( "BURIED_UNSAT_TESTPOSE_PACKSTEPS_1.pdb" ); //DELETE ME.

		utility::vector1< core::conformation::ResidueCOP > residue_vect( pose.total_residue() );

		{
			TR << "TEST 1" << std::endl;
			TR << "---- -" << std::endl;
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				residue_vect[i] = pose.conformation().residue_cop(i);
			}
			core::Real const score1( scoreterm.calculate_energy(residue_vect, 1) );
			TR << "Initial score = " << score1 << std::endl;
			TS_ASSERT_DELTA( score1, 2025, 1e-6 ); // 9 unsatisfieds ((9*5)^2 = 2025)
			core::scoring::EnergyMap totals;
			core::pose::Pose pose1b( pose );
			scoreterm2.finalize_total_energy(pose1b, *sfxn, totals);
			TS_ASSERT_DELTA( score1, totals[core::scoring::buried_unsatisfied_penalty], 1e-6 );
		}

		{
			TR << "TEST 2" << std::endl;
			TR << "---- -" << std::endl;
			core::conformation::ResidueCOP substitution_residue;
			core::conformation::ResidueCOP old_residue( residue_vect[2] );
			for ( core::Size i(1), imax( rotsets->rotamer_set_for_residue(2)->num_rotamers() ); i<=imax; ++i ) {
				substitution_residue = rotsets->rotamer_set_for_residue(2)->rotamer(i);
				if ( substitution_residue->aa() == core::chemical::aa_asp ) break; //Get an asp rotamer.  This turns a donor-acceptor pair into an acceptor-acceptor pair (i.e. removes 1 hbond and creates 1 unsatisfied donor and 1 unsatisfied acceptor).
			}
			residue_vect[2] = substitution_residue;
			core::pose::Pose pose2( pose );
			pose2.replace_residue( 2, *substitution_residue, false );
			//pose2.dump_pdb( "BURIED_UNSAT_TESTPOSE_PACKSTEPS_2.pdb" ); //DELETE ME.
			core::Real const score2( scoreterm.calculate_energy(residue_vect, 2) );
			TR << "Second score = " << score2 << std::endl;
			TS_ASSERT_DELTA( score2, 3025, 1e-6 ); // 11 unsatisfieds ((11*5)^2 = 3025)
			core::scoring::EnergyMap totals;
			scoreterm2.finalize_total_energy(pose2, *sfxn, totals);
			TS_ASSERT_DELTA( score2, totals[core::scoring::buried_unsatisfied_penalty], 1e-6 );
			residue_vect[2] = old_residue; //Reject this one
		}

		core::pose::Pose pose3( pose );
		{
			TR << "TEST 3" << std::endl;
			TR << "---- -" << std::endl;
			core::conformation::ResidueCOP substitution_residue;
			for ( core::Size i(1), imax( rotsets->rotamer_set_for_residue(5)->num_rotamers() ); i<=imax; ++i ) {
				substitution_residue = rotsets->rotamer_set_for_residue(5)->rotamer(i);
				if ( substitution_residue->aa() == core::chemical::aa_asp ) break; //Get an asp rotamer.  This turns a donor-acceptor pair into an acceptor-acceptor pair (i.e. removes 1 hbond and creates 1 unsatisfied donor and 1 unsatisfied acceptor).
			}
			residue_vect[5] = substitution_residue;
			pose3.replace_residue( 5, *substitution_residue, false );
			//pose3.dump_pdb( "BURIED_UNSAT_TESTPOSE_PACKSTEPS_3.pdb" ); //DELETE ME.
			core::Real const score3( scoreterm.calculate_energy(residue_vect, 5) );
			TR << "Third score = " << score3 << std::endl;
			TS_ASSERT_DELTA( score3, 3025, 1e-6 ); // 11 unsatisfieds ((11*5)^2 = 3025)
			core::scoring::EnergyMap totals;
			scoreterm2.finalize_total_energy(pose3, *sfxn, totals);
			TS_ASSERT_DELTA( score3, totals[core::scoring::buried_unsatisfied_penalty], 1e-6 );
			scoreterm.commit_considered_substitution();
			TR << "Accepting substitution 3" << std::endl;
		}

		core::pose::Pose pose4( pose3 );
		{
			TR << "TEST 4" << std::endl;
			TR << "---- -" << std::endl;
			core::conformation::ResidueCOP substitution_residue;
			for ( core::Size i(1), imax( rotsets->rotamer_set_for_residue(8)->num_rotamers() ); i<=imax; ++i ) {
				substitution_residue = rotsets->rotamer_set_for_residue(8)->rotamer(i);
				if ( substitution_residue->aa() == core::chemical::aa_asp ) break; //Get an asp rotamer.  This turns a donor-acceptor pair into an acceptor-acceptor pair (i.e. removes 1 hbond and creates 1 unsatisfied donor and 1 unsatisfied acceptor).
			}
			residue_vect[8] = substitution_residue;
			pose4.replace_residue( 8, *substitution_residue, false );
			//pose4.dump_pdb( "BURIED_UNSAT_TESTPOSE_PACKSTEPS_4.pdb" ); //DELETE ME.
			core::Real const score4( scoreterm.calculate_energy(residue_vect, 8) );
			TR << "Fourth score = " << score4 << std::endl;
			TS_ASSERT_DELTA( score4, 4225, 1e-6 ); // 13 unsatisfieds ((13*5)^2 = 4225)
			core::scoring::EnergyMap totals;
			scoreterm2.finalize_total_energy(pose4, *sfxn, totals);
			TS_ASSERT_DELTA( score4, totals[core::scoring::buried_unsatisfied_penalty], 1e-6 );
			scoreterm.commit_considered_substitution();
			TR << "Accepting substitution 4" << std::endl;
		}

		core::pose::Pose pose5( pose4 );
		{
			TR << "TEST 5" << std::endl;
			TR << "---- -" << std::endl;
			core::conformation::ResidueCOP substitution_residue( pose.conformation().residue_cop(8) );
			core::conformation::ResidueCOP old_residue( residue_vect[8] );
			residue_vect[8] = substitution_residue;
			pose5.replace_residue( 8, *substitution_residue, false );
			//pose5.dump_pdb( "BURIED_UNSAT_TESTPOSE_PACKSTEPS_5.pdb" ); //DELETE ME.
			core::Real const score5( scoreterm.calculate_energy(residue_vect, 8) );
			TR << "Fifth score = " << score5 << std::endl;
			TS_ASSERT_DELTA( score5, 3025, 1e-6 ); // 11 unsatisfieds ((11*5)^2 = 3025)
			core::scoring::EnergyMap totals;
			scoreterm2.finalize_total_energy(pose5, *sfxn, totals);
			TS_ASSERT_DELTA( score5, totals[core::scoring::buried_unsatisfied_penalty], 1e-6 );
			TR << "Rejecting substitution 5" << std::endl;
			residue_vect[8] = old_residue;
		}

		core::pose::Pose pose6( pose4 );
		{
			TR << "TEST 6" << std::endl;
			TR << "---- -" << std::endl;
			core::conformation::ResidueCOP substitution_residue( pose.conformation().residue_cop(5) );
			residue_vect[5] = substitution_residue;
			pose6.replace_residue( 5, *substitution_residue, false );
			//pose6.dump_pdb( "BURIED_UNSAT_TESTPOSE_PACKSTEPS_6.pdb" ); //DELETE ME.
			core::Real const score6( scoreterm.calculate_energy(residue_vect, 5) );
			TR << "Fifth score = " << score6 << std::endl;
			TS_ASSERT_DELTA( score6, 3025, 1e-6 ); // 11 unsatisfieds ((11*5)^2 = 3025)
			core::scoring::EnergyMap totals;
			scoreterm2.finalize_total_energy(pose6, *sfxn, totals);
			TS_ASSERT_DELTA( score6, totals[core::scoring::buried_unsatisfied_penalty], 1e-6 );
			scoreterm.commit_considered_substitution();
			TR << "Accepting substitution 6" << std::endl;
		}

		core::pose::Pose pose7( pose6 );
		{
			TR << "TEST 7" << std::endl;
			TR << "---- -" << std::endl;
			core::conformation::ResidueCOP substitution_residue( pose.conformation().residue_cop(8) );
			residue_vect[8] = substitution_residue;
			pose7.replace_residue( 8, *substitution_residue, false );
			//pose7.dump_pdb( "BURIED_UNSAT_TESTPOSE_PACKSTEPS_7.pdb" ); //DELETE ME.
			core::Real const score7( scoreterm.calculate_energy(residue_vect, 8) );
			TR << "Fifth score = " << score7 << std::endl;
			TS_ASSERT_DELTA( score7, 2025, 1e-6 ); // 9 unsatisfieds ((9*5)^2 = 2025)
			core::scoring::EnergyMap totals;
			scoreterm2.finalize_total_energy(pose7, *sfxn, totals);
			TS_ASSERT_DELTA( score7, totals[core::scoring::buried_unsatisfied_penalty], 1e-6 );
			scoreterm.commit_considered_substitution();
			TR << "Accepting substitution 7" << std::endl;
		}
	}

	/// @brief Test acutally packing a simple pose.
	void test_packing() {
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/guidance_scoreterms/buried_unsat_penalty/graph/hbnet_testcase2.pdb"  );

		core::scoring::methods::EnergyMethodOptions options;
		options.buried_unsatisfied_penalty_burial_threshold(1.0);
		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::buried_unsatisfied_penalty, 1.0 );
		//sfxn->set_weight( core::scoring::hbond_sc, 1.0 );
		sfxn->set_energy_method_options(options);

		//Mutate key residues to asp.
		protocols::simple_moves::MutateResidue mutres2( 2, "ASP" );
		protocols::simple_moves::MutateResidue mutres5( 5, "ASP" );
		protocols::simple_moves::MutateResidue mutres8( 8, "ASP" );
		mutres2.apply(pose);
		mutres5.apply(pose);
		mutres8.apply(pose);

		core::pack::task::PackerTaskOP task( new core::pack::task::PackerTask_(pose));

		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(1))).prevent_repacking();
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(3))).prevent_repacking();
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(4))).prevent_repacking();
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(6))).prevent_repacking();
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(7))).prevent_repacking();
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(9))).prevent_repacking();

		utility::vector1< bool > allowed_aas( core::chemical::num_canonical_aas, false );
		allowed_aas[ core::chemical::aa_asn ] = true;
		allowed_aas[ core::chemical::aa_asp ] = true;
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(2))).restrict_absent_canonical_aas( allowed_aas );
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(5))).restrict_absent_canonical_aas( allowed_aas );
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(8))).restrict_absent_canonical_aas( allowed_aas );

		core::pack::pack_rotamers( pose, *sfxn, task );

		TS_ASSERT_EQUALS( pose.residue_type(2).aa(), core::chemical::aa_asn );
		TS_ASSERT_EQUALS( pose.residue_type(5).aa(), core::chemical::aa_asn );
		TS_ASSERT_EQUALS( pose.residue_type(8).aa(), core::chemical::aa_asn );

		(utility::pointer::static_pointer_cast< core::pack::guidance_scoreterms::buried_unsat_penalty::graph::BuriedUnsatPenaltyGraphContainer const >( pose.energies().data().get_ptr( core::scoring::EnergiesCacheableDataType::BURIED_UNSAT_HBOND_GRAPH ) ))->graph()->provide_pymol_commands_to_show_groups( TR, pose );
		TR.flush();

		//pose.dump_pdb( "BURIED_UNSAT_TESTPOSE_PACKRESULT.pdb" ); //DELETE ME.

	}



};
