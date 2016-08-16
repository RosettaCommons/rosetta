// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/TaskOperations.cxxtest.hh
/// @brief  Tests for TaskOperations reside here.  Happy testing.
/// @author Ben Stranges

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

//Unit headers
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/ReplicateTask.hh>

// project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

// utility headers

//Auto Headers
#include <utility/vector1.hh>


static basic::Tracer TR("test.core.pack.task.operation.TaskOperationsTests");


// --------------- Test Class --------------- //

class TaskOperationsTests : public CxxTest::TestSuite {

	core::pose::Pose pose;

public:

	// --------------- Fixtures --------------- //

	//ctor sets up the pose once
	TaskOperationsTests(){
	}
	virtual ~TaskOperationsTests() {}
	static TaskOperationsTests* createSuite() {
		return new TaskOperationsTests();
	}
	static void destroySuite( TaskOperationsTests *suite ) { delete suite; }

	void setUp() {
		//core_init();
		core_init_with_additional_options( "-mute core.init core.pack.task core.conformation" );
		core::import_pose::pose_from_file( pose, "core/pack/task/resfile_test.pdb" , core::import_pose::PDB_file);
	}
	void tearDown() {}

	// ------------- Helper Functions ------------- //


	// --------------- Test Cases --------------- //

	/// @brief this test applies the DisallowIfNonnative TaskOperation to a PackerTask and verifies that the task
	/// reflects the proper designable residues

	void test_DisallowIfNonnativeTaskOperation() {

		using namespace core;
		using namespace core::pack::task;
		TR << "Running test_DisallowIfNonnativeTaskOperation..." << std::endl;

		core::pack::task::TaskFactory tf;
		operation::DisallowIfNonnativeOP disallow_op( new operation::DisallowIfNonnative() );

		//define options for disallow_op
		std::string noAGP ("AGP");  //don't allow Ala, Gly, Pro
		disallow_op->disallow_aas(noAGP);
		tf.push_back( disallow_op );

		core::pack::task::PackerTaskOP task = tf.create_task_and_apply_taskoperations( pose );
		//TR<< "\n" << *(task) << "\n"<< std::endl; //uncomment to reprint the task

		//now run the test
		test::UTracer UT_disallow("core/pack/task/operation/DisallowIfNonnativeTaskOperation.u");
		UT_disallow << *(task) << std::endl;

	} //end test_DisallowIfNonnative

	/// @brief this test applies the ReplicateTask TaskOperation
	void test_ReplicateTaskTaskOperation(){
		using namespace core;
		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		TR << "Running test_ReplicateTaskTaskOperation..." << std::endl;

		TaskFactoryOP tf( new TaskFactory );
		//define some other operations to do stuffs
		operation::RestrictResidueToRepackingOP nodes_op( new operation::RestrictResidueToRepacking );
		nodes_op->include_residue(1);
		nodes_op->include_residue(2);
		nodes_op->include_residue(4);
		nodes_op->include_residue(31);
		operation::PreventRepackingOP prevent_op( new operation::PreventRepacking );
		prevent_op->include_residue(2);
		prevent_op->include_residue(3);
		prevent_op->include_residue(7);
		//now include one that won't be picked up
		tf->push_back( nodes_op );
		tf->push_back( prevent_op );
		//make the native task
		PackerTaskOP native_task = tf->create_task_and_apply_taskoperations( pose );
		//TR << "NATIVE TASK:" << *(native_task) << std::endl;
		//make a mutation pose that will have a different sequence from the test pose
		//this test makes sure the overall logic is the same even though the sequence is not
		pose::Pose mut_pose;
		TaskFactoryOP mut_tf( new TaskFactory );
		core::import_pose::pose_from_file( mut_pose, "core/pack/task/resfile_test_mut.pdb" , core::import_pose::PDB_file);
		PackerTaskOP mut_task = mut_tf->create_task_and_apply_taskoperations( mut_pose );
		//TR << "MUT TASK BEFORE ReplicateTask: " << *(mut_task)<< std::endl;
		//now make the same logic as the native
		TaskOperationOP mimic_nat_task_op( new operation::ReplicateTask(pose, tf) );
		mimic_nat_task_op->apply( pose, *(mut_task) );
		//TR << "MUT TASK AFTER ReplicateTask: " << *(mut_task)<< std::endl;

		//output for making the UTracers
		//TR<< *(native_task) << *(mut_task) <<std::endl; //uncomment to reprint the task

		//now run the test
		test::UTracer UT_replicate("core/pack/task/operation/ReplicateTaskOperation.u");
		UT_replicate << *(native_task) << *(mut_task) << std::endl;

	} //end of test_ReplicateTaskTaskOperation

}; //end of class deffinition
