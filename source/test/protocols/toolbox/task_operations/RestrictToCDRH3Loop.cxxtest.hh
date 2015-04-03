// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/RestrictNonSurfaceToRepackingOperation.cxxtest.hh
/// @brief  test for the task operation RestrictNonSurfaceToRepackingOperation
/// @author Ron Jacak ronj@unc.edu

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <protocols/toolbox/task_operations/RestrictToCDRH3Loop.hh>

// project headers
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

// utility headers
#include <basic/Tracer.hh>

// C++ headers


static basic::Tracer TR("test.protocols.toolbox.task_operations.RestrictToCDRH3LoopTest");


// --------------- Test Class --------------- //

class RestrictToCDRH3LoopTest : public CxxTest::TestSuite {

	core::pose::Pose pose;

public:

	// --------------- Fixtures --------------- //

	//ctor sets up the pose once
	RestrictToCDRH3LoopTest(){
        core_init();
		core::import_pose::pose_from_pdb( pose, "protocols/util/aHo_antibody.pdb" );
	}
	virtual ~RestrictToCDRH3LoopTest() {}
	static RestrictToCDRH3LoopTest* createSuite() {
		return new RestrictToCDRH3LoopTest();
	}
	static void destroySuite( RestrictToCDRH3LoopTest *suite ) { delete suite; }

	void setUp() {}
	void tearDown() {}

	// ------------- Helper Functions ------------- //
    void perform_relevant_assertions( core::pack::task::PackerTaskOP const & task, core::Size residue_number )
    {
        if ( residue_number >= pose.pdb_info()->pdb2pose( 'H', 107 ) && residue_number <= pose.pdb_info()->pdb2pose( 'H', 138 ) )
        {
            TS_ASSERT_EQUALS( task->residue_task(residue_number).being_packed(), true );
        }
        else
        {
            TS_ASSERT_EQUALS( task->residue_task(residue_number).being_packed(), false );
        }
        TS_ASSERT_EQUALS( task->residue_task(residue_number).being_designed(), false );
    }

	// --------------- Test Cases --------------- //

	/// @brief this test applies the RestrictNonSurfaceToRepackingOperation task operation to a PackerTask and verifies
	/// that all of the residues that have more neighbors than the cutoff have been set to repacking only and that all
	/// the positions with neighbors fewer than the cutoff are designable.

	void test_RestrictToCDRH3Loop() {
    
		
		using protocols::toolbox::task_operations::RestrictToCDRH3Loop;
		using core::pack::task::operation::TaskOperationCOP;
		//using namespace ObjexxFCL::format;

		TR << "Running test_RestrictNonSurfaceToRepackingOperation..." << std::endl;

		core::pack::task::TaskFactory tf;
		tf.push_back( TaskOperationCOP( new RestrictToCDRH3Loop() ) );

		// the following call should work even if the pose hasn't been scored
		core::pack::task::PackerTaskOP task = tf.create_task_and_apply_taskoperations( pose );
        
        for ( core::Size residue_number = 1; residue_number <= pose.total_residue(); ++residue_number )
        {
            perform_relevant_assertions( task, residue_number);
        }
        

	}

};

