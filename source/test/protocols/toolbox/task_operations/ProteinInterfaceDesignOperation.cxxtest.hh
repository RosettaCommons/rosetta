// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/ProteinInterfaceDesignOperation.cxxtest.hh
/// @brief  test for ProteinInterfaceDesignOperation
/// @author Ben Stranges (stranges@unc.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

// Unit header

#include <protocols/toolbox/task_operations/ProteinInterfaceDesignOperation.hh>

// project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>


// C++ headers
//#include <set>
static basic::Tracer TR("test.protocols.toolbox.task_operations.ProteinInterfaceDesignOperationTests");


// --------------- Test Class --------------- //

class ProteinInterfaceDesignOperationTests : public CxxTest::TestSuite {

	core::pose::PoseOP pose;
	core::scoring::ScoreFunctionOP scfxn;

public:

	// --------------- Fixtures --------------- //


	void setUp() {
		core_init();
		//reuse for comparison with Interface class
		pose = core::pose::PoseOP( new core::pose::Pose() );
		scfxn = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() );
		core::import_pose::pose_from_pdb( *pose, "core/conformation/dock_in.pdb" );

		//need to score the pose to find the interface in this case
		scfxn->set_weight( core::scoring::fa_atr, 0.8  );
		scfxn->set_weight( core::scoring::fa_rep, 0.44 );
		(*scfxn)(*pose);

	}

	void tearDown() {}

	// ------------- Helper Functions ------------- //


	// --------------- Test Cases --------------- //

	//test for default operation settings
	void test_ProteinInterfaceDesignOperation() {

		using namespace core::pack::task;
		using protocols::toolbox::task_operations::ProteinInterfaceDesignOperation;
		using core::pack::task::operation::TaskOperationCOP;
		TaskFactory PID_factory;
		PID_factory.push_back( TaskOperationCOP( new ProteinInterfaceDesignOperation() ) );
		TR << "Running test_ProteinInterfaceDesignOperation" << std::endl;
		test::UTracer UT_PID("protocols/toolbox/task_operations/ProteinInterfaceDesignOperation.u");
		UT_PID << *(PID_factory.create_task_and_apply_taskoperations( *pose ) ) << std::endl;
	} // end test_ProteinInterfaceDesignOperation


	//test for designing all aas at interface
	void test_ProteinInterfaceDesignOperation_all() {

		using namespace core::pack::task;
		using namespace protocols::toolbox::task_operations;
		ProteinInterfaceDesignOperationOP pid_op( new ProteinInterfaceDesignOperation() );
		//allows design of Cys, Gly, Pro at all positions
		pid_op->allow_all_aas( true );
		pid_op->design_all_aas( true );
		TaskFactory PID_factory;
		PID_factory.push_back( pid_op );
		TR << "Running test_ProteinInterfaceDesignOperation_all" << std::endl;
		test::UTracer UT_PID("protocols/toolbox/task_operations/ProteinInterfaceDesignOperation_all.u");
		UT_PID << *(PID_factory.create_task_and_apply_taskoperations( *pose ) ) << std::endl;
	} // end test_ProteinInterfaceDesignOperation_all


};
