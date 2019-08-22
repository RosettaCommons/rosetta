// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/SequenceSymmetricAnnealer.cxxtest.hh
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/pack/annealer/SequenceSymmetricAnnealer.hh>

#include <core/chemical/AA.hh>

#include <utility/graph/Graph.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/import_pose/import_pose.hh>

// Test headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <core/pack/task/operation/KeepSequenceSymmetry.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/task_operations/SetIGTypeOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/simple_moves/VirtualRootMover.hh>

#include <list>
#include <basic/Tracer.hh>

static basic::Tracer TR( "core.pack.annealer.SequenceSymmetricAnnealer.cxxtest" );

using namespace core;
using namespace core::pack;

class SequenceSymmetricAnnealerTests : public CxxTest::TestSuite {
public:

	void setUp() {
		core_init();
	}

	//Pack a 3-chain protein, check that all chains end up with the same sequence
	void test_chain_similarity_after_packing() {
		TR << "starting test_chain_similarity_after_packing" << std::endl;
		std::list< core::pack::task::operation::TaskOperationCOP > operations;
		operations.emplace_back( utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >() );

		{//for speedup
			utility::vector1< std::string > base_types = { "ALA", "GLY", "SER" };
			operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types ) );
		}

		protocols::minimization_packing::PackRotamersMover packer;
		packer.initialize_task_factory_with_operations( operations );

		pose::PoseOP pose =
			import_pose::pose_from_file( "core/scoring/symmetry/2akf.pdb" );
		TS_ASSERT_EQUALS( pose->num_chains(), 3 );
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 2 ) );
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 3 ) );

		auto const old_sequence = pose->chain_sequence( 1 );
		packer.apply( * pose );
		TS_ASSERT( pose->chain_sequence( 1 ) != old_sequence );//Failure here is a failure of the unit test, not of the protocol
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 2 ) );//Failer here is a failure of the protocol
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 3 ) );//Failer here is a failure of the protocol

		TR << "New sequences: " << pose->chain_sequence( 1 ) << " " << pose->chain_sequence( 2 ) << " " << pose->chain_sequence( 3 ) << std::endl;
	}

	//Pack a 3-chain protein, check that all chains end up with the same sequence
	// This one uses the commandline interface
	void test_chain_similarity_after_packing_cli() {
		TR << "starting test_chain_similarity_after_packing_cli" << std::endl;
		//-sequence_symmetric_annealer 1
		core_init_with_additional_options( "-sequence_symmetric_annealer 1" );
		std::list< core::pack::task::operation::TaskOperationCOP > operations;
		operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::InitializeFromCommandline >() );

		{//for speedup
			utility::vector1< std::string > base_types = { "ALA", "GLY", "SER" };
			operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types ) );
		}

		protocols::minimization_packing::PackRotamersMover packer;
		packer.initialize_task_factory_with_operations( operations );

		pose::PoseOP pose =
			import_pose::pose_from_file( "core/scoring/symmetry/2akf.pdb" );
		TS_ASSERT_EQUALS( pose->num_chains(), 3 );
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 2 ) );
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 3 ) );

		auto const old_sequence = pose->chain_sequence( 1 );
		packer.apply( * pose );
		TS_ASSERT( pose->chain_sequence( 1 ) != old_sequence );//Failure here is a failure of the unit test, not of the protocol
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 2 ) );//Failer here is a failure of the protocol
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 3 ) );//Failer here is a failure of the protocol

		TR << "New sequences: " << pose->chain_sequence( 1 ) << " " << pose->chain_sequence( 2 ) << " " << pose->chain_sequence( 3 ) << std::endl;
	}

	//Pack a 3-chain protein, check that all chains end up with the same sequence
	// This one has a virtual root
	void test_chain_similarity_with_virtual_root() {
		TR << "starting test_chain_similarity_with_virtual_root" << std::endl;
		//-sequence_symmetric_annealer 1
		core_init_with_additional_options( "-sequence_symmetric_annealer 1" );
		std::list< core::pack::task::operation::TaskOperationCOP > operations;
		operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::InitializeFromCommandline >() );

		{//for speedup
			utility::vector1< std::string > base_types = { "ALA", "GLY", "SER" };
			operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types ) );
		}

		protocols::minimization_packing::PackRotamersMover packer;
		packer.initialize_task_factory_with_operations( operations );

		pose::PoseOP pose =
			import_pose::pose_from_file( "core/scoring/symmetry/2akf.pdb" );
		TS_ASSERT_EQUALS( pose->num_chains(), 3 );
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 2 ) );
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 3 ) );

		protocols::simple_moves::VirtualRootMover vrm;
		vrm.apply( * pose );
		//TS_ASSERT_EQUALS( pose->num_chains(), 4 ); //Not true, adds "X" residue to last chain

		auto const old_sequence = pose->chain_sequence( 1 );
		packer.apply( * pose );
		TS_ASSERT( pose->chain_sequence( 1 ) != old_sequence );//Failure here is a failure of the unit test, not of the protocol
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 2 ) );//Failer here is a failure of the protocol

		//Add an X for the virtual root
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ) + "X", pose->chain_sequence( 3 ) );//Failer here is a failure of the protocol

		TR << "New sequences: " << pose->chain_sequence( 1 ) << " " << pose->chain_sequence( 2 ) << " " << pose->chain_sequence( 3 ) << std::endl;
	}


	//Pack a 3-chain protein, check that all chains end up with the same sequence
	//In this version, the first chain STARTS with a different sequence than the other two
	void test_recover_symmetry() {
		TR << "starting test_recover_symmetry" << std::endl;
		std::list< core::pack::task::operation::TaskOperationCOP > operations;
		operations.emplace_back( utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >() );

		{//for speedup
			utility::vector1< std::string > base_types = { "ALA", "GLY", "SER" };
			operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types ) );
		}

		pose::PoseOP pose =
			import_pose::pose_from_file( "core/scoring/symmetry/2akf.pdb" );

		std::string const seq_for_chain1 =
			"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
		///"SAGASAGASAGASAGASAGASAGASAGASAGA";

		auto op = utility::pointer::make_shared< core::pack::task::operation::RestrictResidueToRepacking >( );
		for ( core::Size resid = 1; resid <= pose->chain_sequence( 1 ).size(); ++resid ) {
			protocols::simple_moves::MutateResidue mut( resid, seq_for_chain1[ resid - 1 ] );
			mut.apply( * pose );
			op->include_residue( resid );
		}
		operations.push_back( op );

		TS_ASSERT_EQUALS( seq_for_chain1, pose->chain_sequence( 1 ) );
		TS_ASSERT( seq_for_chain1 != pose->chain_sequence( 2 ) );
		TS_ASSERT( seq_for_chain1 != pose->chain_sequence( 3 ) );

		protocols::minimization_packing::PackRotamersMover packer;
		packer.initialize_task_factory_with_operations( operations );
		packer.apply( * pose );

		TS_ASSERT_EQUALS( seq_for_chain1, pose->chain_sequence( 1 ) );//Failer here is a failure of the protocol
		TS_ASSERT_EQUALS( seq_for_chain1, pose->chain_sequence( 2 ) );//Failer here is a failure of the protocol
		TS_ASSERT_EQUALS( seq_for_chain1, pose->chain_sequence( 3 ) );//Failer here is a failure of the protocol

		TR << "New sequences: " << pose->chain_sequence( 1 ) << " " << pose->chain_sequence( 2 ) << " " << pose->chain_sequence( 3 ) << std::endl;
	}

	//Pack a 3-chain protein, check that all chains end up with the same sequence
	//In this version, chain 1 is missing its first residue and chain3 is missing its last residue
	void test_unequal_chain_lengths(){
		TR << "starting test_unequal_chain_lengths" << std::endl;
		std::list< core::pack::task::operation::TaskOperationCOP > operations;
		operations.emplace_back( utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >() );

		{//for speedup
			utility::vector1< std::string > base_types = { "ALA", "VAL", "LYS" };
			operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types ) );
		}

		protocols::minimization_packing::PackRotamersMover packer;
		packer.initialize_task_factory_with_operations( operations );

		std::string const expected_chain2_seq = "VSRLEEDVRNLNAIVQKLQERLDRLEETVQAK";

		pose::PoseOP pose =
			import_pose::pose_from_file( "core/pack/annealer/2akf.clipped.pdb" );
		TS_ASSERT_EQUALS( pose->num_chains(), 3 );
		TS_ASSERT_EQUALS( pose->chain_sequence( 2 ), expected_chain2_seq );
		TS_ASSERT_EQUALS( pose->chain_sequence( 2 ), "V" + pose->chain_sequence( 1 ) );
		TS_ASSERT_EQUALS( pose->chain_sequence( 2 ), pose->chain_sequence( 3 ) + "K" );

		auto const old_sequence = pose->chain_sequence( 1 );
		packer.apply( * pose );
		TS_ASSERT( pose->chain_sequence( 1 ) != old_sequence );//Failure here is a failure of the unit test, not of the protocol

		TS_ASSERT_EQUALS( pose->chain_sequence( 2 ), "V" + pose->chain_sequence( 1 ) );
		TS_ASSERT_EQUALS( pose->chain_sequence( 2 ), pose->chain_sequence( 3 ) + "K" );

		TR << "New sequences: _" << pose->chain_sequence( 1 ) << " " << pose->chain_sequence( 2 ) << " " << pose->chain_sequence( 3 ) << "_" << std::endl;
	}

	//Provide the packer an impossible problem
	//This is a branch off of test_unequal_chain_lengths
	//This problem is impossible because the packer has repack-only positions (because some chains are missing residues at that position) but we are forcing mutations to alanine
	//COMMENTING OUT BECAUSE I CAN'T HANDLE ERRORS ELEGANTLY
	/*
	void NEGATIVE_test_fail_if_impossible(){
	std::list< core::pack::task::operation::TaskOperationCOP > operations;
	operations.emplace_back( utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >() );

	utility::vector1< std::string > base_types = { "ALA" };
	operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types ) );

	protocols::minimization_packing::PackRotamersMover packer;
	packer.initialize_task_factory_with_operations( operations );

	std::string const expected_chain2_seq = "VSRLEEDVRNLNAIVQKLQERLDRLEETVQAK";

	pose::PoseOP pose =
	import_pose::pose_from_file( "core/pack/annealer/2akf.clipped.pdb" );
	TS_ASSERT_EQUALS( pose->num_chains(), 3 );
	TS_ASSERT_EQUALS( pose->chain_sequence( 2 ), expected_chain2_seq );
	TS_ASSERT_EQUALS( pose->chain_sequence( 2 ), "V" + pose->chain_sequence( 1 ) );
	TS_ASSERT_EQUALS( pose->chain_sequence( 2 ), pose->chain_sequence( 3 ) + "K" );

	bool error = false;
	try {
	packer.apply( * pose );
	} catch( ... ){
	error = true;
	}
	TS_ASSERT( error );
	}
	*/

};
