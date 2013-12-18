// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/devel/zn_hash/ZnHash.cxxtest.hh
/// @brief  test suite for devel::zn_hash::ZnHash
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Unit headers
#include <devel/mmt_msd/MMTMinPackingJob.hh>

/// Project headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/scmin/SidechainStateAssignment.hh>

// --------------- Test Class --------------- //

class MMTMinPackingJobTests : public CxxTest::TestSuite {

public:

	typedef core::Size   Size;
	typedef core::Real   Real;



	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	void
	initialize_minpacking_job( devel::mmt_msd::MMTMinPackingJobOP job ) {
		core::pose::Pose pose = create_trpcage_ideal_pose();
		core::scoring::ScoreFunctionOP sfxn = core::scoring::getScoreFunction();

		core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( ii != 11 ) {
				task->nonconst_residue_task( ii ).prevent_repacking();
			} else {
				utility::vector1< bool > restrict_absent( core::chemical::num_canonical_aas, false );
				restrict_absent[ core::chemical::aa_phe ] = true;
				task->nonconst_residue_task( ii ).restrict_absent_canonical_aas( restrict_absent );
			}
		}

		job->set_pose( pose );
		job->set_sfxn( *sfxn );
		job->set_packer_task( *task );

	}

	devel::mmt_msd::MMTMinPackingJobOP
	create_and_initialize_minpacking_job() {
		devel::mmt_msd::MMTMinPackingJobOP job = new devel::mmt_msd::MMTMinPackingJob;
		initialize_minpacking_job( job );
		return job;
	}

	/// @details Make sure when you set the pose, sfxn, and task, that the desired properties
	/// for the class are actually set and reflect their settings.
	void test_initialize_minpacking_job() {
		devel::mmt_msd::MMTMinPackingJobOP job =  create_and_initialize_minpacking_job();
		core::scoring::ScoreFunctionOP dft_sfxn = core::scoring::getScoreFunction();

		TS_ASSERT( job->has_pose() );
		TS_ASSERT( job->has_sfxn() );
		TS_ASSERT( job->has_task() );

		// if any of the above three TS_ASSERT statement fails, we have a major problem
		if ( ! job->has_pose() || ! job->has_sfxn() || ! job->has_task() ) return;

		TS_ASSERT( job->get_pose().total_residue() == 20 );
		TS_ASSERT( job->get_sfxn().weights() == dft_sfxn->weights() );
		TS_ASSERT( job->get_task()->num_to_be_packed() == 1 );
		TS_ASSERT( job->get_task()->design_any() );
		TS_ASSERT( job->get_task()->being_packed( 11 ) );
		TS_ASSERT( job->get_task()->being_designed( 11 ) );
		TS_ASSERT( ! job->best_assignment_exists() );
	}

	void test_run_quick_minpacking() {
		devel::mmt_msd::MMTMinPackingJobOP job = create_and_initialize_minpacking_job();
		core::pose::Pose trpcage = create_trpcage_ideal_pose();
		TS_ASSERT( trpcage.residue( 11 ).aa() == core::chemical::aa_gly );

		job->go();

		TS_ASSERT( job->optimization_complete() );
		TS_ASSERT( job->running_time() > 0 );

		TS_ASSERT( job->best_assignment_exists() );
		TS_ASSERT( ! job->get_best_assignment().any_unassigned() );
		TS_ASSERT( ! job->has_pose() );
		TS_ASSERT( ! job->has_sfxn() );
		TS_ASSERT( ! job->has_task() );

		initialize_minpacking_job( job );

		TS_ASSERT( job->has_pose() );
		TS_ASSERT( job->has_sfxn() );
		TS_ASSERT( job->has_task() );

		// now go and rebuild the data needed to reconstruct the final pose
		job->setup();
		job->update_pose( trpcage );
		TS_ASSERT( trpcage.residue( 11 ).aa() == core::chemical::aa_phe );
	}

};
