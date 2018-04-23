// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pack/guidance_scoreterms/buried_unsat_penalty/BuriedUnsatPenaltySymmetricTests.cxxtest.hh
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
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/ResidueLevelTask.hh>
#include <core/pack/task/ResidueLevelTask_.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/symmetry/SymmetricEnergies.hh>

// Protocols Headers, for convenience
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/symmetry/SetupForSymmetryMover.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("BuriedUnsatPenaltySymmetricTests");


class BuriedUnsatPenaltySymmetricTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}

	/// @brief Simple test: just score a pose.
	/// @details This pose was created with Scott's HBNet, and should have fully satisfied hydrogen bonds.
	void test_simple_symmetric_scoring(){
		core::pose::Pose asymmpose;
		core::import_pose::pose_from_file( asymmpose, "core/pack/guidance_scoreterms/buried_unsat_penalty/graph/hbnet_testcase1.pdb"  );
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/guidance_scoreterms/buried_unsat_penalty/graph/hbnet_testcase1_symm.pdb"  );
		{ //Set up symmetry:
			protocols::symmetry::SetupForSymmetryMover setupsymm;
			setupsymm.process_symmdef_file( "core/pack/guidance_scoreterms/buried_unsat_penalty/C3.symm" );
			setupsymm.apply(pose);
		}

		core::scoring::ScoreFunction asymmsfxn;
		asymmsfxn.set_weight( core::scoring::buried_unsatisfied_penalty, 1.0 );

		try {
			asymmsfxn(asymmpose);
		} catch( utility::excn::Exception & exception ) {
			TR << "Exception caught in unit test: " << exception.msg() << std::endl;
			TS_ASSERT(false);
		} catch ( ... ) {
			TR << "Unhandled exception caught in unit test!" << std::endl;
			TS_ASSERT(false);
		}

		core::scoring::symmetry::SymmetricScoreFunction sfxn;
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

		core::Real const EXPECTED_ENERGY( std::pow( 9*5, 2 ) ); //Nine unsats.  Could also #define this.  Meh.  It's a unit test.
		TS_ASSERT_DELTA( pose.energies().total_energy(), EXPECTED_ENERGY, 1e-6 );
		TS_ASSERT_DELTA( asymmpose.energies().total_energy(), EXPECTED_ENERGY, 1e-6 );
		TS_ASSERT_DELTA( pose.energies().total_energy(), asymmpose.energies().total_energy(), 1e-6 );

		pose.dump_pdb("BURIED_UNSAT_PENALTY_SYMMETRIC_SIMPLE_SCORING.pdb"); //DELETE ME
		{
			//The following nonsense is just to dump out the groups:
			core::scoring::methods::EnergyMethodOptions options;
			core::pack::guidance_scoreterms::buried_unsat_penalty::BuriedUnsatPenalty scoreterm( options );
			core::pack::task::PackerTaskOP task( new core::pack::task::PackerTask_(pose) );
			for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
				(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(i))).restrict_to_repacking();
			}
			task->request_symmetrize_by_intersection();
			core::pack::task::PackerTaskOP symmtask( core::pack::make_new_symmetric_PackerTask_by_requested_method( pose, task ) );

			core::pack::rotamer_set::symmetry::SymmetricRotamerSetsOP rotsets( utility::pointer::make_shared< core::pack::rotamer_set::symmetry::SymmetricRotamerSets >() );
			rotsets->set_task( symmtask );
			rotsets->initialize_pose_for_rotsets_creation(pose);
			utility::graph::GraphOP packer_neighbor_graph( core::pack::create_packer_graph( pose, sfxn, symmtask ) );
			rotsets->build_rotamers(pose, sfxn, packer_neighbor_graph);
			rotsets->prepare_sets_for_packing(pose, sfxn);
			scoreterm.set_up_residuearrayannealableenergy_for_packing( pose, *rotsets, sfxn );
			scoreterm.provide_pymol_commands_to_show_groups( TR, pose );
			TR.flush();
		}

		protocols::simple_moves::MutateResidue mutres( 54, "ASP" ); //Mutate a hydrogen bond-forming residue.
		mutres.apply(pose);
		// This should result in 6 acceptors unsatisfied.
		sfxn(pose); //Rescore
		TS_ASSERT_DELTA( pose.energies().total_energy(), static_cast< core::Real >( std::pow( 15*5, 2 ) ), 1e-6 );
	}

	/// @brief Test penalty graph setup for packing.
	void test_symmetric_setup_for_packing() {
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/guidance_scoreterms/buried_unsat_penalty/graph/hbnet_testcase1_symm.pdb"  );
		{ //Set up symmetry:
			protocols::symmetry::SetupForSymmetryMover setupsymm;
			setupsymm.process_symmdef_file( "core/pack/guidance_scoreterms/buried_unsat_penalty/C3.symm" );
			setupsymm.apply(pose);
		}

		core::scoring::methods::EnergyMethodOptions options;
		core::pack::guidance_scoreterms::buried_unsat_penalty::BuriedUnsatPenalty scoreterm( options );
		core::scoring::symmetry::SymmetricScoreFunctionOP sfxn( new core::scoring::symmetry::SymmetricScoreFunction );
		sfxn->set_weight(core::scoring::buried_unsatisfied_penalty, 1.0);
		//sfxn->set_weight( core::scoring::hbond_sc, 1.0 );

		(*sfxn)(pose);

		core::pack::task::PackerTaskOP task( new core::pack::task::PackerTask_(pose) );
		for ( core::Size i(1), imax(pose.total_residue()); i<=imax; ++i ) {
			(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(i))).restrict_to_repacking();
		}

		task->request_symmetrize_by_intersection();
		core::pack::task::PackerTaskOP symmtask( core::pack::make_new_symmetric_PackerTask_by_requested_method( pose, task ) );

		core::pack::rotamer_set::symmetry::SymmetricRotamerSetsOP rotsets( new core::pack::rotamer_set::symmetry::SymmetricRotamerSets );
		rotsets->set_task( symmtask );
		rotsets->initialize_pose_for_rotsets_creation(pose);
		utility::graph::GraphOP packer_neighbor_graph( core::pack::create_packer_graph( pose, *sfxn, symmtask ) );
		rotsets->build_rotamers(pose, *sfxn, packer_neighbor_graph);
		rotsets->prepare_sets_for_packing(pose, *sfxn);

		scoreterm.set_up_residuearrayannealableenergy_for_packing( pose, *rotsets, *sfxn );
		core::pack::guidance_scoreterms::buried_unsat_penalty::graph::BuriedUnsatPenaltyGraph const & bunsat_graph( *(scoreterm.unsat_graph_) );
		TS_ASSERT_EQUALS( rotsets->total_residue(), pose.total_residue() ); //Should be true
		for ( core::Size i(1), imax(75 /*Pose is 75 residues per asymm unit*/); i<=imax; ++i ) {
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

	/// @brief Test acutally packing a simple pose.
	void test_symmetric_packing() {
		core::pose::Pose pose;
		core::import_pose::pose_from_file( pose, "core/pack/guidance_scoreterms/buried_unsat_penalty/graph/hbnet_testcase2_symm.pdb"  );
		{ //Set up symmetry:
			protocols::symmetry::SetupForSymmetryMover setupsymm;
			setupsymm.process_symmdef_file( "core/pack/guidance_scoreterms/buried_unsat_penalty/C3.symm" );
			setupsymm.apply(pose);
		}

		core::scoring::methods::EnergyMethodOptions options;
		options.buried_unsatisfied_penalty_burial_threshold(1.0);
		core::scoring::symmetry::SymmetricScoreFunctionOP sfxn( new core::scoring::symmetry::SymmetricScoreFunction );
		sfxn->set_weight( core::scoring::buried_unsatisfied_penalty, 1.0 );
		//sfxn->set_weight( core::scoring::hbond_sc, 1.0 );
		sfxn->set_energy_method_options(options);

		//Mutate key residues to asp.
		protocols::simple_moves::MutateResidue mutres2( 2, "ASP" );
		mutres2.apply(pose);

		core::pack::task::PackerTaskOP task( new core::pack::task::PackerTask_(pose));

		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(1))).restrict_to_repacking();
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(3))).restrict_to_repacking();

		utility::vector1< bool > allowed_aas( core::chemical::num_canonical_aas, false );
		allowed_aas[ core::chemical::aa_asn ] = true;
		allowed_aas[ core::chemical::aa_asp ] = true;
		(static_cast< core::pack::task::ResidueLevelTask_ & >(task->nonconst_residue_task(2))).restrict_absent_canonical_aas( allowed_aas );

		task->request_symmetrize_by_intersection();
		core::pack::task::PackerTaskOP symmtask( core::pack::make_new_symmetric_PackerTask_by_requested_method( pose, task ) );

		core::pack::pack_rotamers( pose, *sfxn, symmtask );

		TS_ASSERT_EQUALS( pose.residue_type(2).aa(), core::chemical::aa_asn );
		TS_ASSERT_EQUALS( pose.residue_type(5).aa(), core::chemical::aa_asn );
		TS_ASSERT_EQUALS( pose.residue_type(8).aa(), core::chemical::aa_asn );

		(utility::pointer::static_pointer_cast< core::pack::guidance_scoreterms::buried_unsat_penalty::graph::BuriedUnsatPenaltyGraphContainer const >( pose.energies().data().get_ptr( core::scoring::EnergiesCacheableDataType::BURIED_UNSAT_HBOND_GRAPH ) ))->graph()->provide_pymol_commands_to_show_groups( TR, pose );
		TR.flush();

		pose.dump_pdb( "BURIED_UNSAT_TESTPOSE_SYMM_PACKRESULT.pdb" ); //DELETE ME.

	}



};
