// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResfileReader.cxxtest.hh
/// @brief  test suite for resfile reader
/// @author Steven Lewis

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <core/types.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/RotamerSampleOptions.hh>

#include <core/pack/task/ResfileReader.hh>

#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh> // for a handful of options that should be read
#include <basic/options/keys/in.OptionKeys.gen.hh> // for a handful of options that should not be read
#include <basic/options/option.hh>
#include <basic/options/util.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/backtrace.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKey.hh>

// C++ headers
#include <string>
#include <sstream>



// --------------- Test Class --------------- //

class PackerTaskTests : public CxxTest::TestSuite {

private:
	core::pose::PoseOP pose;
	core::pack::task::PackerTaskOP task;

public:

	typedef utility::keys::VariantKey< utility::options::OptionKey > VariantOptionKey;


	// --------------- Fixtures --------------- //

	void setUp() {
		core_init();
		pose = create_trpcage_ideal_poseop();
		task = core::pack::task::TaskFactory::create_packer_task( *pose );
	}

	void tearDown() {
		pose.reset();
		task.reset();
	}


	void test_packer_task_default_behavior_design_all_positions() {
		TS_ASSERT( pose->total_residue() == task->total_residue() );
		if ( pose->total_residue() != task->total_residue() ) return;
		TS_ASSERT( task->num_to_be_packed() == task->total_residue() );
		TS_ASSERT( task->design_any() );
		for ( core::Size ii = 1; ii <= pose->total_residue(); ++ii ) {
			TS_ASSERT( task->pack_residue( ii ) );
			TS_ASSERT( task->being_packed( ii ) );
			TS_ASSERT( task->design_residue( ii ) );
			TS_ASSERT( task->being_designed( ii ) );
		}
	}

	/// make sure that after a single call to restrict_to_repacking (and nothing else since construction)
	/// that every residue in the task is set to repacking, but that none of them are set to being
	/// designed.
	void test_packer_task_restrict_to_repacking() {
		task->restrict_to_repacking();

		TS_ASSERT( pose->total_residue() == task->total_residue() );
		if ( pose->total_residue() != task->total_residue() ) return;
		TS_ASSERT( task->num_to_be_packed() == task->total_residue() );
		TS_ASSERT( ! task->design_any() );
		for ( core::Size ii = 1; ii <= pose->total_residue(); ++ii ) {
			TS_ASSERT( task->pack_residue( ii ) );
			TS_ASSERT( task->being_packed( ii ) );
			TS_ASSERT( ! task->design_residue( ii ) );
			TS_ASSERT( ! task->being_designed( ii ) );
		}
	}

	/// @brief Test the logic for decing which residues to allow based on the original residue type: namely
	/// that the "name3" criterion is used to decide whether two residue type objects reflect the same underlying
	/// amino acid
	void test_packer_task_name3_decision_making() {
		pose = create_test_in_pdb_poseop(); // need a PDB that has a histadine so we can test for the presence of both HIS and HIS_D
		task = core::pack::task::TaskFactory::create_packer_task( *pose );

		task->restrict_to_repacking();

		TS_ASSERT( pose->total_residue() == task->total_residue() );
		if ( pose->total_residue() != task->total_residue() ) return;

		for ( core::Size ii = 1; ii <= pose->total_residue(); ++ii ) {
			/// the list of residue types should not be empty
			TS_ASSERT( task->residue_task( ii ).allowed_residue_types_begin() != task->residue_task( ii ).allowed_residue_types_end() );
			for ( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter
					iter = task->residue_task( ii ).allowed_residue_types_begin(),
					iter_end = task->residue_task( ii ).allowed_residue_types_end();
					iter != iter_end; ++iter ) {
				TS_ASSERT( (*iter)->name3() == pose->residue_type(ii).name3() );
			}
		}

	}

	/// @brief Test the logic inside the fix_his_tautomer method.
	void test_packer_task_fix_his_tautomer() {
		pose = create_test_in_pdb_poseop(); // need a PDB that has a histadine
		task = core::pack::task::TaskFactory::create_packer_task( *pose );

		task->restrict_to_repacking();

		TS_ASSERT( pose->residue_type( 6 ).aa() == core::chemical::aa_his );
		TS_ASSERT( task->residue_task( 6 ).allowed_residue_types().size() == 2 );

		// 1. test that fix_his_tautomer disables the HIS_D type, since the
		// original test_in.pdb pose is protonated on NE2
		task->nonconst_residue_task( 6 ).or_fix_his_tautomer( true );
		TS_ASSERT( task->residue_task( 6 ).allowed_residue_types().size() == 1 );
		TS_ASSERT( (* task->residue_task( 6 ).allowed_residue_types_begin() )->name() == "HIS" );

		// 2. test that this operation is commutative, in that, it obeys an "or" like behavior;
		// setting fix_his_tautomer to false after first setting it to true doesn't re-enable
		// HIS_D
		task->nonconst_residue_task( 6 ).or_fix_his_tautomer( false );
		TS_ASSERT( task->residue_task( 6 ).allowed_residue_types().size() == 1 );
		TS_ASSERT( (* task->residue_task( 6 ).allowed_residue_types_begin() )->name() == "HIS" );

	}

	void test_packer_task_or_ex1_commutativity() {
		/// should start out "false"
		TS_ASSERT( ! task->residue_task( 1 ).ex1() );

		task->nonconst_residue_task(1).or_ex1( true );
		TS_ASSERT( task->residue_task( 1 ).ex1() );

		// now try to set ex1 to false and observe that it remains true.
		task->nonconst_residue_task(1).or_ex1( false );
		TS_ASSERT( task->residue_task( 1 ).ex1() );
	}

	void test_packer_task_or_ex2_commutativity() {
		/// should start out "false"
		TS_ASSERT( ! task->residue_task( 1 ).ex2() );

		task->nonconst_residue_task(1).or_ex2( true );
		TS_ASSERT( task->residue_task( 1 ).ex2() );

		// now try to set ex2 to false and observe that it remains true.
		task->nonconst_residue_task(1).or_ex2( false );
		TS_ASSERT( task->residue_task( 1 ).ex2() );
	}

	void test_packer_task_or_ex3_commutativity() {
		/// should start out "false"
		TS_ASSERT( ! task->residue_task( 1 ).ex3() );

		task->nonconst_residue_task(1).or_ex3( true );
		TS_ASSERT( task->residue_task( 1 ).ex3() );

		// now try to set ex3 to false and observe that it remains true.
		task->nonconst_residue_task(1).or_ex3( false );
		TS_ASSERT( task->residue_task( 1 ).ex3() );
	}

	void test_packer_task_or_ex4_commutativity() {
		/// should start out "false"
		TS_ASSERT( ! task->residue_task( 1 ).ex4() );

		task->nonconst_residue_task(1).or_ex4( true );
		TS_ASSERT( task->residue_task( 1 ).ex4() );

		// now try to set ex4 to false and observe that it remains true.
		task->nonconst_residue_task(1).or_ex4( false );
		TS_ASSERT( task->residue_task( 1 ).ex4() );
	}

	void test_packer_task_or_ex1aro_commutativity() {
		/// should start out "false"
		TS_ASSERT( ! task->residue_task( 1 ).ex1aro() );

		task->nonconst_residue_task(1).or_ex1aro( true );
		TS_ASSERT( task->residue_task( 1 ).ex1aro() );

		// now try to set ex1aro to false and observe that it remains true.
		task->nonconst_residue_task(1).or_ex1aro( false );
		TS_ASSERT( task->residue_task( 1 ).ex1aro() );
	}

	void test_packer_task_or_ex2aro_commutativity() {
		/// should start out "false"
		TS_ASSERT( ! task->residue_task( 1 ).ex2aro() );

		task->nonconst_residue_task(1).or_ex2aro( true );
		TS_ASSERT( task->residue_task( 1 ).ex2aro() );

		// now try to set ex1aro to false and observe that it remains true.
		task->nonconst_residue_task(1).or_ex2aro( false );
		TS_ASSERT( task->residue_task( 1 ).ex2aro() );
	}

	void test_packer_task_and_extrachi_cutoff() {
		TS_ASSERT( (int) task->residue_task(1).extrachi_cutoff() == core::pack::task::EXTRACHI_CUTOFF_LIMIT );

		/// now try to decrease the extrachi cutoff, and it will successfully be decreased.
		task->nonconst_residue_task(1).and_extrachi_cutoff( 10 );
		TS_ASSERT( task->residue_task(1).extrachi_cutoff() == 10 );

		/// now try to increase the extrachi_cutoff, and it will not increase
		task->nonconst_residue_task(1).and_extrachi_cutoff( 12 );
		TS_ASSERT( task->residue_task(1).extrachi_cutoff() == 10 );
	}

	void test_packer_task_update_commutative_restrict_to_repacking1() {
		core::pack::task::PackerTaskOP task2 = core::pack::task::TaskFactory::create_packer_task( *pose );
		task2->restrict_to_repacking();
		task->update_commutative( *task2 );

		TS_ASSERT( pose->total_residue() == task->total_residue() );
		if ( pose->total_residue() != task->total_residue() ) return;
		TS_ASSERT( task->num_to_be_packed() == task->total_residue() );
		TS_ASSERT( ! task->design_any() );
		for ( core::Size ii = 1; ii <= pose->total_residue(); ++ii ) {
			TS_ASSERT( task->pack_residue( ii ) );
			TS_ASSERT( task->being_packed( ii ) );
			TS_ASSERT( ! task->design_residue( ii ) );
			TS_ASSERT( ! task->being_designed( ii ) );
		}
	}

	void test_packer_task_update_commutative_restrict_to_repacking2() {
		core::pack::task::PackerTaskOP task2 = core::pack::task::TaskFactory::create_packer_task( *pose );
		task->restrict_to_repacking();
		task->update_commutative( *task2 );

		TS_ASSERT( pose->total_residue() == task->total_residue() );
		if ( pose->total_residue() != task->total_residue() ) return;
		TS_ASSERT( task->num_to_be_packed() == task->total_residue() );
		TS_ASSERT( ! task->design_any() );
		for ( core::Size ii = 1; ii <= pose->total_residue(); ++ii ) {
			TS_ASSERT( task->pack_residue( ii ) );
			TS_ASSERT( task->being_packed( ii ) );
			TS_ASSERT( ! task->design_residue( ii ) );
			TS_ASSERT( ! task->being_designed( ii ) );
		}
	}

	void test_PackerTask_initialize_from_options_and_list_options_read_are_in_sync()
	{

		utility::options::OptionKeyList packer_task_options;
		core::pack::task::PackerTask::list_options_read( packer_task_options );
		// cursory check of some of the options known to be read by the ImportPoseOptions
		TS_ASSERT_EQUALS( std::count( packer_task_options.begin(), packer_task_options.end(), VariantOptionKey( basic::options::OptionKeys::packing::linmem_ig )), 1 );
		TS_ASSERT_EQUALS( std::count( packer_task_options.begin(), packer_task_options.end(), VariantOptionKey( basic::options::OptionKeys::packing::fix_his_tautomer )), 1 );
		TS_ASSERT_EQUALS( std::count( packer_task_options.begin(), packer_task_options.end(), VariantOptionKey( basic::options::OptionKeys::packing::max_rotbump_energy )), 1 );

		// cursory check that not all options are in here, because that would be weird
		TS_ASSERT_EQUALS( std::count( packer_task_options.begin(), packer_task_options.end(), VariantOptionKey( basic::options::OptionKeys::in::file::s )), 0 );

		utility::options::OptionCollectionCOP packer_task_option_collection =
			basic::options::read_subset_of_global_option_collection( packer_task_options );

		// Now, try to create an ImportPoseOptions from the local option collection
		try {
			set_throw_on_next_assertion_failure(); // just in case
			task->initialize_from_options( *packer_task_option_collection );
		} catch ( ... ) {
			TS_ASSERT( false ); // we screwed the pooch
		}
	}

	void test_PackerTask_init_from_options_actually_reads_option_collection()
	{

		utility::options::OptionKeyList packer_task_options;
		core::pack::task::PackerTask::list_options_read( packer_task_options );

		// Now drop one of the options from the list, one which always gets read, and make sure that
		// when we call the initialize_from_options function, that an assertion failure occurs
		packer_task_options.remove( VariantOptionKey( basic::options::OptionKeys::packing::fix_his_tautomer ));

		utility::options::OptionCollectionCOP packer_task_option_collection =
			basic::options::read_subset_of_global_option_collection( packer_task_options );

		// Now, try to create an ImportPoseOptions from the local option collection
		try {
			set_throw_on_next_assertion_failure();
			task->initialize_from_options( *packer_task_option_collection );
			TS_ASSERT( false ); // we screwed the pooch
		} catch ( ... ) {
			// good -- if we don't list an option that we're going to read, then
			// an exception will be thrown / an assertion failure will get triggered
			TS_ASSERT( true );
		}

	}

	void test_ResidueLevelTask_initialize_from_options_and_list_options_read_are_in_sync()
	{

		utility::options::OptionKeyList res_lvl_task_options;
		core::pack::task::ResidueLevelTask::list_options_read( res_lvl_task_options );
		// cursory check of some of the options known to be read by the ImportPoseOptions
		TS_ASSERT_EQUALS( std::count( res_lvl_task_options.begin(), res_lvl_task_options.end(), VariantOptionKey( basic::options::OptionKeys::packing::ex1::ex1 )), 1 );
		TS_ASSERT_EQUALS( std::count( res_lvl_task_options.begin(), res_lvl_task_options.end(), VariantOptionKey( basic::options::OptionKeys::packing::use_input_sc )), 1 );
		TS_ASSERT_EQUALS( std::count( res_lvl_task_options.begin(), res_lvl_task_options.end(), VariantOptionKey( basic::options::OptionKeys::packing::extrachi_cutoff )), 1 );

		// cursory check that not all options are in here, because that would be weird
		TS_ASSERT_EQUALS( std::count( res_lvl_task_options.begin(), res_lvl_task_options.end(), VariantOptionKey( basic::options::OptionKeys::in::file::s )), 0 );

		utility::options::OptionCollectionCOP res_lvl_task_option_collection =
			basic::options::read_subset_of_global_option_collection( res_lvl_task_options );

		// Now, try to create an ImportPoseOptions from the local option collection
		try {
			set_throw_on_next_assertion_failure(); // just in case
			task->nonconst_residue_task( 5 ).initialize_from_options( *res_lvl_task_option_collection );
		} catch ( ... ) {
			TS_ASSERT( false ); // we screwed the pooch
		}
	}

	void test_ResidueLevelTask_init_from_options_actually_reads_option_collection()
	{

		utility::options::OptionKeyList res_lvl_task_options;
		core::pack::task::PackerTask::list_options_read( res_lvl_task_options );

		// Now drop one of the options from the list, one which always gets read, and make sure that
		// when we call the initialize_from_options function, that an assertion failure occurs
		res_lvl_task_options.remove( VariantOptionKey( basic::options::OptionKeys::packing::ex1::ex1 ));

		utility::options::OptionCollectionCOP res_lvl_task_option_collection =
			basic::options::read_subset_of_global_option_collection( res_lvl_task_options );

		// Now, try to create an ImportPoseOptions from the local option collection
		try {
			set_throw_on_next_assertion_failure();
			task->nonconst_residue_task( 5 ).initialize_from_options( *res_lvl_task_option_collection );
			TS_ASSERT( false ); // we screwed the pooch
		} catch ( ... ) {
			// good -- if we don't list an option that we're going to read, then
			// an exception will be thrown / an assertion failure will get triggered
			TS_ASSERT( true );
		}

	}

};//end class
