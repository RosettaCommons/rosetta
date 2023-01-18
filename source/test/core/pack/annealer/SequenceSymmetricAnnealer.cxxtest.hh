// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/annealer/SequenceSymmetricAnnealer.cxxtest.hh
/// @author Jack Maguire
/// @author Updated by Tim Neary, timdot10@gmail.com

// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/pack/annealer/SequenceSymmetricAnnealer.fwd.hh>
#include <core/pack/palette/CustomBaseTypePackerPalette.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


#include <core/import_pose/import_pose.hh>

// Test headers
#include <test/core/init_util.hh>

//Auto Headers
#include <utility/vector1.hh>

#include <core/pack/task/operation/KeepSequenceSymmetry.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/task_operations/SetIGTypeOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/simple_moves/VirtualRootMover.hh>
#include <protocols/symmetry/SetupForSequenceSymmetryMover.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/pack/task/xml_util.hh> // For TASK_OPERATIONS_TAG
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>

#include <list>
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh> // AUTO IWYU For Pose

static basic::Tracer TR( "core.pack.annealer.SequenceSymmetricAnnealer.cxxtest" );

using namespace core;
using namespace core::pack;

class SequenceSymmetricAnnealerTests : public CxxTest::TestSuite {
private:

	/// @brief Convenience function to aid in clarity for tests
	/// Builds task operations used to speed up computation of IG
	std::list< core::pack::task::operation::TaskOperationCOP >
	setup_speedup_task_ops() {
		std::list< core::pack::task::operation::TaskOperationCOP > operations;

		{//for speedup
			utility::vector1< std::string > base_types = { "ALA", "GLY", "SER"};
			operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types ) );
			operations.emplace_back( utility::pointer::make_shared< protocols::task_operations::SetIGTypeOperation >( false, false, false, true ) );
		}
		return operations;
	}

	std::list< core::pack::task::operation::TaskOperationCOP >
	setup_speedup__ncaa_task_ops() {
		std::list< core::pack::task::operation::TaskOperationCOP > operations;

		{//for speedup
			utility::vector1< std::string > base_types = { "NVL", "DNVL" };
			operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types ) );
			operations.emplace_back( utility::pointer::make_shared< protocols::task_operations::SetIGTypeOperation >( false, false, false, true ) );
		}
		return operations;
	}

	/// @brief Convenience function to import a test pose and determine if it is correct
	core::pose::PoseOP
	import_and_test_pose() {
		pose::PoseOP pose =
			import_pose::pose_from_file( "core/scoring/symmetry/2akf.pdb" );
		TS_ASSERT_EQUALS( pose->num_chains(), 3 );
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 2 ) );
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 3 ) );
		return pose;
	}

public:

	void setUp() {
		core_init();
	}

	//Pack a 3-chain protein, check that all chains end up with the same sequence
	void test_chain_similarity_after_packing() {
		TR << "starting test_chain_similarity_after_packing" << std::endl;
		auto operations = setup_speedup_task_ops();
		auto kss_top = utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >();
		kss_top->set_prefix_name( "test" );
		operations.push_back( kss_top );

		auto ss_mop = utility::pointer::make_shared< protocols::symmetry::SetupForSequenceSymmetryMover >();
		ss_mop->set_prefix_name( kss_top->prefix_name() );

		// Add chain residue selectors
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "1" ) );
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "2" ) );
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "3" ) );

		auto pose = import_and_test_pose();

		ss_mop->apply( * pose );
		protocols::minimization_packing::PackRotamersMover packer;
		packer.initialize_task_factory_with_operations( operations );

		auto const old_sequence = pose->chain_sequence( 1 );
		packer.apply( * pose );
		TS_ASSERT( pose->chain_sequence( 1 ) != old_sequence );//Failure here is a failure of the unit test, not of the protocol
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 2 ) );//Failer here is a failure of the protocol
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 3 ) );//Failer here is a failure of the protocol

		// TR << "New sequences: " << pose->chain_sequence( 1 ) << " " << pose->chain_sequence( 2 ) << " " << pose->chain_sequence( 3 ) << std::endl;
	}

	//Pack a 3-chain protein, only allow noncanonical amino acids at one position, check that all chains end up with the same sequence
	void test_chain_similarity_after_packing_ncaa() {
		TR << "starting test_chain_similarity_after_packing_ncaa" << std::endl;
		auto operations = setup_speedup__ncaa_task_ops();
		auto kss_top = utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >();
		kss_top->set_prefix_name( "test" );
		operations.push_back( kss_top );
		//utility::vector1< std::string > ncaa_base_types = { "NVL", "DNVL" };
		//operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( ncaa_base_types ) );

		auto ss_mop = utility::pointer::make_shared< protocols::symmetry::SetupForSequenceSymmetryMover >();
		ss_mop->set_prefix_name( kss_top->prefix_name() );

		// Add chain residue selectors
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "1" ) );
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "2" ) );
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "3" ) );

		// Add PackerPalette with L- and D-Norvaline
		using namespace core::pack::palette;
		core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() );
		core::pack::task::TaskFactoryOP task = utility::pointer::make_shared< core::pack::task::TaskFactory >();
		CustomBaseTypePackerPaletteOP ncaa_pp = utility::pointer::make_shared< CustomBaseTypePackerPalette >();
		ncaa_pp->parse_additional_residue_types( "NVL,DNVL" );
		task->set_packer_palette( ncaa_pp );
		for ( const auto & op : operations ) {
			task->push_back(op);
		};
		auto pose = import_and_test_pose();
		auto const old_sequence = pose->annotated_sequence( 1 );

		ss_mop->apply( * pose );
		protocols::minimization_packing::PackRotamersMover packer(sfxn, task->create_task_and_apply_taskoperations( * pose));
		packer.apply( * pose);

		TS_ASSERT( pose->annotated_sequence( 1 ) != old_sequence );//Failure here is a failure of the unit test, not of the protocol
		TS_ASSERT_EQUALS( pose->annotated_sequence( 1 ), pose->annotated_sequence( 2 ) );//Failer here is a failure of the protocol
		TS_ASSERT_EQUALS( pose->annotated_sequence( 1 ), pose->annotated_sequence( 3 ) );//Failer here is a failure of the protocol

		TR << "New sequences: " << pose->annotated_sequence( 1 ) << " " << pose->annotated_sequence( 2 ) << " " << pose->annotated_sequence( 3 ) << std::endl;
	}

	// Probably more trouble than its worth to implement a cli method for this...
	//Pack a 3-chain protein, check that all chains end up with the same sequence
	// This one uses the commandline interface
	/*
	void NEGATIVE_test_chain_similarity_after_packing_cli() {
	TR << "starting test_chain_similarity_after_packing_cli" << std::endl;
	//-sequence_symmetric_annealer 1
	core_init_with_additional_options( "-sequence_symmetric_annealer 1" );
	std::list< core::pack::task::operation::TaskOperationCOP > operations;
	operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::InitializeFromCommandline >() );
	operations.emplace_back( utility::pointer::make_shared< protocols::task_operations::SetIGTypeOperation >( false, false, false, true ) );

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
	*/

	//Pack a 3-chain protein, check that all chains end up with the same sequence
	// This one has a virtual root
	void test_chain_similarity_with_virtual_root() {
		TR << "starting test_chain_similarity_with_virtual_root" << std::endl;
		auto operations = setup_speedup_task_ops();
		auto kss_top = utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >();
		kss_top->set_prefix_name( "test" );
		operations.push_back( kss_top );

		auto ss_mop = utility::pointer::make_shared< protocols::symmetry::SetupForSequenceSymmetryMover >();
		ss_mop->set_prefix_name( kss_top->prefix_name() );

		// Add chain residue selectors
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "1" ) );
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "2" ) );
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "3" ) );

		auto pose = import_and_test_pose();

		ss_mop->apply( * pose );
		protocols::minimization_packing::PackRotamersMover packer;
		packer.initialize_task_factory_with_operations( operations );

		protocols::simple_moves::VirtualRootMover vrm;
		vrm.apply( * pose );
		TS_ASSERT_EQUALS( pose->num_chains(), 4 );

		auto const old_sequence = pose->chain_sequence( 1 );
		packer.apply( * pose );
		TS_ASSERT_DIFFERS( pose->chain_sequence( 1 ), old_sequence );  // Failure here is a failure of the unit test, not of the protocol
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 2 ) );  // Failure here is a failure of the unit test, not of the protocol
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 3 ) );  // Failure here is a failure of the unit test, not of the protocol
		TS_ASSERT_EQUALS( pose->chain_sequence( 4 ), "X" );

		//Add an X for the virtual root
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 3 ) );//Failer here is a failure of the protocol

		// TR << "New sequences: " << pose->chain_sequence( 1 ) << " " << pose->chain_sequence( 2 ) << " " << pose->chain_sequence( 3 ) << std::endl;
	}

	//Pack a 3-chain protein, check that all chains end up with the same sequence
	// This test will test to see whether linked residues linked using different regions give the correct result
	// i.e. linking chain 1 and 2, followed by linking chain 2 and 3 is functionally identical to linking 1,2,3 in one go.
	void test_chain_similiarity_with_alternative_logic() {
		TR << "starting test_chain_similiarity_with_alternative_logic" << std::endl;
		auto operations = setup_speedup_task_ops();
		auto kss_top = utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >();
		kss_top->set_prefix_name( "test" );
		operations.push_back( kss_top );

		auto ss_mop = utility::pointer::make_shared< protocols::symmetry::SetupForSequenceSymmetryMover >();
		ss_mop->set_prefix_name( kss_top->prefix_name() );

		// Add chain residue selectors
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "1" ) );
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "2" ) );
		ss_mop->add_residue_selector( 1, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "2" ) );
		ss_mop->add_residue_selector( 1, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "3" ) );

		auto pose = import_and_test_pose();

		ss_mop->apply( * pose );
		protocols::minimization_packing::PackRotamersMover packer;
		packer.initialize_task_factory_with_operations( operations );

		auto const old_sequence = pose->chain_sequence( 1 );
		packer.apply( * pose );
		TS_ASSERT( pose->chain_sequence( 1 ) != old_sequence );// Failure here is a failure of the unit test, not of the protocol
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 2 ) );// Failure here is a failure of the protocol
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 3 ) );// Failure here is a failure of the protocol

		// TR << "New sequences: " << pose->chain_sequence( 1 ) << " " << pose->chain_sequence( 2 ) << " " << pose->chain_sequence( 3 ) << std::endl;
	}

	//Pack a 3-chain protein, check that all chains end up with the same sequence
	//In this version, the first chain STARTS with a different sequence than the other two
	void test_recover_symmetry() {
		TR << "starting test_recover_symmetry" << std::endl;
		auto operations = setup_speedup_task_ops();
		auto kss_top = utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >();
		kss_top->set_prefix_name( "test" );
		operations.push_back( kss_top );

		auto ss_mop = utility::pointer::make_shared< protocols::symmetry::SetupForSequenceSymmetryMover >();
		ss_mop->set_prefix_name( kss_top->prefix_name() );

		// Add chain residue selectors
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "1" ) );
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "2" ) );
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "3" ) );

		auto pose = import_and_test_pose();

		std::string const seq_for_chain1( pose->chain_sequence( 1 ).size(), 'A' );

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

		ss_mop->apply( * pose );
		protocols::minimization_packing::PackRotamersMover packer;
		packer.initialize_task_factory_with_operations( operations );
		packer.apply( * pose );

		TS_ASSERT_EQUALS( seq_for_chain1, pose->chain_sequence( 1 ) );//Failer here is a failure of the protocol
		TS_ASSERT_EQUALS( seq_for_chain1, pose->chain_sequence( 2 ) );//Failer here is a failure of the protocol
		TS_ASSERT_EQUALS( seq_for_chain1, pose->chain_sequence( 3 ) );//Failer here is a failure of the protocol

		// TR << "New sequences: " << pose->chain_sequence( 1 ) << " " << pose->chain_sequence( 2 ) << " " << pose->chain_sequence( 3 ) << std::endl;
	}

	// From the new functionaility, this will now not work as originally intended -> all but the last two entries will be identical
	// If this type of method is attempted from scripts rosetta will exit before packing as selectors are not equal.
	//Pack a 3-chain protein, check that all chains end up with the same sequence
	//In this version, chain 1 is missing its first residue and chain3 is missing its last residue
	/*
	void NEGATIVE_test_unequal_chain_lengths(){
	TR << "starting test_unequal_chain_lengths" << std::endl;
	std::list< core::pack::task::operation::TaskOperationCOP > operations;
	operations.emplace_back( utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >() );
	operations.emplace_back( utility::pointer::make_shared< protocols::task_operations::SetIGTypeOperation >( false, false, false, true ) );

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
	*/

	// Pack using 3 linked residues, where each selector defines different residues
	// Linked residues 1-10 to 10-20 respectively
	void test_index_residue_selectors() {
		TR << "starting test_index_residue_selectors" << std::endl;
		auto operations = setup_speedup_task_ops();
		auto kss_top = utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >();
		kss_top->set_prefix_name( "test" );
		operations.push_back( kss_top );

		auto ss_mop = utility::pointer::make_shared< protocols::symmetry::SetupForSequenceSymmetryMover >();
		ss_mop->set_prefix_name( kss_top->prefix_name() );

		// Add chain residue selectors
		utility::vector1< core::Size > const res_idx_1 = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
		utility::vector1< core::Size > const res_idx_2 = { 11,12,13,14,15,16,17,18,19,20 };
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( res_idx_1 ) );
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( res_idx_2 ) );

		auto pose = import_and_test_pose();

		ss_mop->apply( * pose );
		protocols::minimization_packing::PackRotamersMover packer;
		packer.initialize_task_factory_with_operations( operations );

		auto const old_sequence = pose->chain_sequence( 1 );
		// TR << "PACKING" << std::endl;
		packer.apply( * pose );
		std::string const seq = pose->sequence();
		TS_ASSERT( pose->chain_sequence( 1 ) != old_sequence );//Failure here is a failure of the unit test, not of the protocol
		TS_ASSERT_EQUALS( seq.substr(0, 10), seq.substr(10, 10) ); // Failure here is failure of protocol

		// TR << "Linked residues from 1-10 and 11-20." << std::endl;
		// TR << "New sequences: " << pose->chain_sequence( 1 ) << " " << pose->chain_sequence( 2 ) << " " << pose->chain_sequence( 3 ) << std::endl;
	}

	// Pack using 3 linked residues, where there is overlap in the index selector.
	// Linked residues 1-10 to 3-13 respectively
	// The expected output is to have odd and even residues linked.
	void test_overlapping_index_residue_selectors() {
		TR << "starting test_overlapping_index_residue_selectors" << std::endl;
		auto operations = setup_speedup_task_ops();
		auto kss_top = utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >();
		kss_top->set_prefix_name( "test" );
		operations.push_back( kss_top );

		auto ss_mop = utility::pointer::make_shared< protocols::symmetry::SetupForSequenceSymmetryMover >();
		ss_mop->set_prefix_name( kss_top->prefix_name() );

		// Add chain residue selectors
		utility::vector1< core::Size > const res_idx_1 = { 1,2,3,4,5,6,7, 8, 9,10 };
		utility::vector1< core::Size > const res_idx_2 = { 3,4,5,6,7,8,9,10,11,12 };
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( res_idx_1 ) );
		ss_mop->add_residue_selector( 0, utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >( res_idx_2 ) );

		auto pose = import_and_test_pose();

		ss_mop->apply( * pose );
		protocols::minimization_packing::PackRotamersMover packer;
		packer.initialize_task_factory_with_operations( operations );

		auto const old_sequence = pose->chain_sequence( 1 );
		// TR << "PACKING" << std::endl;
		packer.apply( * pose );
		std::string const seq = pose->sequence();
		TS_ASSERT( pose->chain_sequence( 1 ) != old_sequence );//Failure here is a failure of the unit test, not of the protocol
		for ( core::Size ii = 0; ii < 10; ++ii ) {
			TS_ASSERT_EQUALS( seq[ii], seq[ii + 2] ); // Linked res ii and ii + 2, e.g. 1 and 3
		}

		// TR << "All odd and even numbers for residues from 1 to 12 should be linked." << std::endl;
		// TR << "New sequences: " << pose->chain_sequence( 1 ) << " " << pose->chain_sequence( 2 ) << " " << pose->chain_sequence( 3 ) << std::endl;
	}

	//Provide the packer an impossible problem
	//This is a branch off of test_unequal_chain_lengths
	//This problem is impossible because the packer has repack-only positions (because some chains are missing residues at that position) but we are forcing mutations to alanine
	//COMMENTING OUT BECAUSE I CAN'T HANDLE ERRORS ELEGANTLY
	/*
	void NEGATIVE_test_fail_if_impossible(){
	std::list< core::pack::task::operation::TaskOperationCOP > operations;
	operations.emplace_back( utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >() );
	operations.emplace_back( utility::pointer::make_shared< protocols::task_operations::SetIGTypeOperation >( false, false, false, true ) );

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

	// This test is used to determine whether linked residues with different acceptable residue types are properly handled.
	// Here chain B will only be able to mutate to Ala and chain B can mutate to any of Ala,Gly,Ser.
	// The result should be a poly-Ala chain for both.
	void test_different_accepted_restypes_for_linked_residues() {
		TR << "starting test_different_accepted_restypes_for_linked_residues" << std::endl;
		std::list< core::pack::task::operation::TaskOperationCOP > operations;
		auto kss_top = utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >();
		kss_top->set_prefix_name( "test" );
		operations.emplace_back( kss_top );

		auto chain1 = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "1" );
		auto chain2 = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "2" );
		auto chain3 = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "3" );

		auto ss_mop = utility::pointer::make_shared< protocols::symmetry::SetupForSequenceSymmetryMover >();
		ss_mop->set_prefix_name( kss_top->prefix_name() );
		ss_mop->add_residue_selector( 0, chain1 );
		ss_mop->add_residue_selector( 0, chain2 );

		utility::vector1< std::string > base_types1 = { "ALA", "GLY", "SER" };
		utility::vector1< std::string > base_types2 = { "ALA" };
		operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types1, chain1 ) );
		operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types2, chain2 ) );
		operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types2, chain3 ) );
		operations.emplace_back( utility::pointer::make_shared< protocols::task_operations::SetIGTypeOperation >( false, false, false, true ) );

		auto pose = import_and_test_pose();

		auto const old_sequence = pose->chain_sequence( 1 );
		ss_mop->apply( * pose );
		protocols::minimization_packing::PackRotamersMover packer;
		packer.initialize_task_factory_with_operations( operations );

		packer.apply( * pose );
		TS_ASSERT( pose->chain_sequence( 1 ) != old_sequence );//Failure here is a failure of the unit test, not of the protocol
		std::string const poly_ala( pose->chain_sequence( 1 ).size(), 'A' );
		TS_ASSERT_EQUALS( poly_ala, pose->chain_sequence( 1 ) );
		TS_ASSERT_EQUALS( poly_ala, pose->chain_sequence( 2 ) );
	}

	// This test is for checking the error handling of the Annealer. In this instance an error is expected as
	// no  common residue types are present for any of the linked residues.
	void test_fail_in_packing_from_no_common_res_types() {
		TR << "starting test_different_accepted_restypes_for_linked_residues" << std::endl;
		std::list< core::pack::task::operation::TaskOperationCOP > operations;
		auto kss_top = utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >();
		kss_top->set_prefix_name( "test" );
		operations.emplace_back( kss_top );

		auto chain1 = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "1" );
		auto chain2 = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "2" );
		auto chain3 = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "3" );

		auto ss_mop = utility::pointer::make_shared< protocols::symmetry::SetupForSequenceSymmetryMover >();
		ss_mop->set_prefix_name( kss_top->prefix_name() );
		ss_mop->add_residue_selector( 0, chain1 );
		ss_mop->add_residue_selector( 0, chain2 );

		utility::vector1< std::string > base_types1 = { "GLY", "SER" }; // No common residue types
		utility::vector1< std::string > base_types2 = { "ALA" };  // No common residue types
		operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types1, chain1 ) );
		operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types2, chain2 ) );
		operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types2, chain3 ) ); // For speedup
		operations.emplace_back( utility::pointer::make_shared< protocols::task_operations::SetIGTypeOperation >( false, false, false, true ) );

		auto pose = import_and_test_pose();

		auto const old_sequence = pose->chain_sequence( 1 );
		ss_mop->apply( * pose );
		protocols::minimization_packing::PackRotamersMover packer;
		packer.initialize_task_factory_with_operations( operations );

		try {
			packer.apply( * pose );
			TR << "SequenceSymmetricAnnealer should not have completed without throwing an exception, no common residue types for linked res were present." << std::endl;
			TS_ASSERT( false );
		} catch( utility::excn::Exception & e ) {
			TS_ASSERT( true );
		}
	}

	//Test that NVL and DNVL are not treated as interchangeable (this also tests the more general case, that different NCAAs are not treated as interchangeable)

	void test_fail_in_packing_from_no_common_res_types_L_and_D() {
		TR << "starting test_different_accepted_restypes_L_and_D_for_linked_residues" << std::endl;
		std::list< core::pack::task::operation::TaskOperationCOP > operations;
		auto kss_top = utility::pointer::make_shared< task::operation::KeepSequenceSymmetry >();
		kss_top->set_prefix_name( "test" );
		operations.emplace_back( kss_top );

		auto chain1 = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "1" );
		auto chain2 = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "2" );
		auto chain3 = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "3" );

		auto ss_mop = utility::pointer::make_shared< protocols::symmetry::SetupForSequenceSymmetryMover >();
		ss_mop->set_prefix_name( kss_top->prefix_name() );
		ss_mop->add_residue_selector( 0, chain1 );
		ss_mop->add_residue_selector( 0, chain2 );

		utility::vector1< std::string > base_types1 = { "NVL", "SER" }; // No common residue types
		utility::vector1< std::string > base_types2 = { "DNVL" };  // No common residue types
		operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types1, chain1 ) );
		operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types2, chain2 ) );
		operations.emplace_back( utility::pointer::make_shared< core::pack::task::operation::RestrictToSpecifiedBaseResidueTypes >( base_types2, chain3 ) ); // For speedup
		operations.emplace_back( utility::pointer::make_shared< protocols::task_operations::SetIGTypeOperation >( false, false, false, true ) );

		using namespace core::pack::palette;
		core::scoring::ScoreFunctionOP sfxn( core::scoring::get_score_function() );
		core::pack::task::TaskFactoryOP task = utility::pointer::make_shared< core::pack::task::TaskFactory >();
		CustomBaseTypePackerPaletteOP ncaa_pp = utility::pointer::make_shared< CustomBaseTypePackerPalette >();
		ncaa_pp->parse_additional_residue_types( "NVL,DNVL" );
		task->set_packer_palette( ncaa_pp );
		for ( const auto & op : operations ) {
			task->push_back(op);
		};
		auto pose = import_and_test_pose();

		auto const old_sequence = pose->chain_sequence( 1 );
		ss_mop->apply( * pose );
		protocols::minimization_packing::PackRotamersMover packer(sfxn, task->create_task_and_apply_taskoperations( * pose));

		try {
			packer.apply( * pose );
			TR << "SequenceSymmetricAnnealer should not have completed without throwing an exception, no common residue types for linked res were present." << std::endl;
			TS_ASSERT( false );
		} catch( utility::excn::Exception & e ) {
			TS_ASSERT( true );
		}
	}

	// Tests for validation of a proper xml schema. Ensures that each of the relevant parts are properly used.
	// Uses three tests to ensure that other linked tags must be parsed properly too.
	void test_xml_validation() {
		TR << "starting test_xml_validation" << std::endl;

		std::string kss_xml =
			"<KeepSequenceSymmetry name=\"test_kss\" setting=\"true\"/>";
		std::string sss_xml =
			"<SetupForSequenceSymmetry name=\"test_sss\" sequence_symmetry_behaviour=\"test_kss\">\n"
			"\t<SequenceSymmetry residue_selectors=\"cA,cB\"/>\n"
			"</SetupForSequenceSymmetry>";
		std::string sss_fail_xml =
			"<SetupForSequenceSymmetry name=\"test_sss\" sequence_symmetry_behaviour=\"test_kss\">\n"
			"\t<SequenceSymmetry residue_selectors=\"cA,cI\"/>\n"
			"</SetupForSequenceSymmetry>";

		auto pose = import_and_test_pose();

		auto res1_sel = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "1" );
		auto res2_sel = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "2" );
		auto resI_sel = utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >(
			utility::vector1< core::Size > { 1, 2, 3, 4 } );

		// Create tags to test XML parsing
		utility::tag::TagCOP kss_tag   = utility::tag::Tag::create( kss_xml );
		utility::tag::TagCOP sss_tag   = utility::tag::Tag::create( sss_xml );
		utility::tag::TagCOP sss_tag_f = utility::tag::Tag::create( sss_fail_xml );
		basic::datacache::DataMap dm;

		auto kss_top = utility::pointer::make_shared< core::pack::task::operation::KeepSequenceSymmetry >();
		auto sss_mov = utility::pointer::make_shared< protocols::symmetry::SetupForSequenceSymmetryMover >();

		dm.add( "ResidueSelector", "cA", res1_sel ); // Needed for parsing SSS
		dm.add( "ResidueSelector", "cB", res2_sel ); // Needed for parsing SSS
		dm.add( "ResidueSelector", "cI", resI_sel ); // Needed for parsing SSS

		// Check KeepSequenceSymmetry task op is parsed correctly
		try {
			kss_top->parse_tag( kss_tag, dm );
			TS_ASSERT( true );
		} catch( utility::excn::Exception & e ) {
			TR << "KeepSequenceSymmetry should not have been parsed incorrectly." << std::endl;
			TS_ASSERT( false );
		}
		// As SetupForSequenceSymmetry require KeepSequenceSymmetry this is expected to fail
		try {
			sss_mov->parse_my_tag( sss_tag, dm );
			TR << "SetupForSequenceSymmetryMover should not have been parsed. No KeepSequenceSymmetry was present in the datamap." << std::endl;
			TS_ASSERT( false ); // Should fail
		} catch( utility::excn::Exception & e ) {
			dm.add( core::pack::task::TASK_OPERATIONS_TAG, "test_kss", kss_top );
			TS_ASSERT( true );
		}
		// This should fail due to differently sized residue selectors
		try {
			sss_mov->parse_my_tag( sss_tag_f, dm );
			sss_mov->apply( * pose );
			TR << "SetupForSequenceSymmetryMover should not have parsed. The given residue selectors were of unequal length." << std::endl;
			TS_ASSERT( false );
		} catch( utility::excn::Exception & e ) {
			TS_ASSERT( true ); // Should fail in this instance
		}
		// Now KeepSequenceSymmetry is in the datamap it should pass correctly.
		try {
			sss_mov->parse_my_tag( sss_tag, dm );
			sss_mov->apply( * pose );
			TS_ASSERT( true );
		} catch( utility::excn::Exception & e ) {
			TR << "SetupForSequenceSymmetry was parsed incorrectly, it should have been accepted." << std::endl;
			TS_ASSERT( false ); // This should not fail anymore.
		}
	}

	// Tests for validation of a proper xml schema. It parses an xml and ensures that the resulting protein is packed correctly.
	void test_initialise_annealer_from_xml() {
		TR << "starting test_initialise_annealer_from_xml" << std::endl;

		std::string kss_xml =
			"<KeepSequenceSymmetry name=\"test_kss\" setting=\"true\"/>";
		std::string sss_xml =
			"<SetupForSequenceSymmetry name=\"test_sss\" sequence_symmetry_behaviour=\"test_kss\">\n"
			"\t<SequenceSymmetry residue_selectors=\"cA,cB\"/>\n"
			"</SetupForSequenceSymmetry>";

		auto pose = import_and_test_pose();

		auto res1_sel = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "1" );
		auto res2_sel = utility::pointer::make_shared< core::select::residue_selector::ChainSelector >( "2" );

		// Create tags to test XML parsing
		utility::tag::TagCOP kss_tag = utility::tag::Tag::create( kss_xml );
		utility::tag::TagCOP sss_tag = utility::tag::Tag::create( sss_xml );
		basic::datacache::DataMap dm;

		auto kss_top = utility::pointer::make_shared< core::pack::task::operation::KeepSequenceSymmetry >();
		auto sss_mov = utility::pointer::make_shared< protocols::symmetry::SetupForSequenceSymmetryMover >();

		dm.add( "ResidueSelector", "cA", res1_sel ); // Needed for parsing SSS
		dm.add( "ResidueSelector", "cB", res2_sel ); // Needed for parsing SSS

		kss_top->parse_tag( kss_tag, dm );
		dm.add( core::pack::task::TASK_OPERATIONS_TAG, "test_kss", kss_top );
		sss_mov->parse_my_tag( sss_tag, dm );
		sss_mov->apply( * pose );

		// Now to pack the protein to assess whether the mover has been properly processed.
		auto operations = setup_speedup_task_ops();
		operations.push_back( kss_top );

		protocols::minimization_packing::PackRotamersMover packer;
		packer.initialize_task_factory_with_operations( operations );

		auto const old_sequence = pose->chain_sequence( 1 );
		packer.apply( * pose );
		TS_ASSERT( pose->chain_sequence( 1 ) != old_sequence );// Failure here is a failure of the unit test, not of the protocol
		TS_ASSERT_EQUALS( pose->chain_sequence( 1 ), pose->chain_sequence( 2 ) );//Failer here is a failure of the protocol
	}

};
