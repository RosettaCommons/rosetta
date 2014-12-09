// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/TaskAwareMinMover.cxxtest.hh
/// @brief  test for TaskAwareMinMover
/// @author Steven Lewis

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Utilities
#include <core/types.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/kinematics/MoveMap.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

//Auto Headers
#include <core/id/AtomID_Mask.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


// --------------- Test Class --------------- //

class TaskAwareMinMoverTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //

	void setUp() {
		core_init();
	}

	void tearDown() {
	}

	// ------------- Helper Functions ------------- //


	// --------------- Test Cases --------------- //

	///@details The purpose of TAMinMover is to cause MinMover to respect a PackerTask, and minimize whatever sidechains are packable.  This test sets up three ((TaskFactory, TAMinMover), (MinMover, MoveMap)) quadruplets where the task-aware part should have the same effect as the hard-coded MoveMap part.
	void test_TaskAwareMinMover() {

		core::pose::Pose const backup(create_twores_1ubq_pose());
		core::pose::Pose p1(backup), p2(backup);

		//////////////////////TaskFactory(s)
		using namespace core::pack::task;
		TaskFactoryOP task_factory_all( new TaskFactory() );
		operation::PreventRepackingOP prop_all( new operation::PreventRepacking );
		prop_all->include_residue(1); prop_all->include_residue(2);
		task_factory_all->push_back( prop_all );

		TaskFactoryOP task_factory_one( new TaskFactory() );
		operation::PreventRepackingOP prop_one( new operation::PreventRepacking );
		prop_one->include_residue(1);
		task_factory_one->push_back( prop_one );

		TaskFactoryOP task_factory( new TaskFactory() );
		operation::PreventRepackingOP prop( new operation::PreventRepacking );
		task_factory->push_back( prop );

		///////////////////////starting MoveMaps
		core::kinematics::MoveMapOP mm_start( new core::kinematics::MoveMap() );
		mm_start->set_bb(false); mm_start->set_jump(false); mm_start->set_chi(false);

		/////////////////////////minimizer movers/////////////////////////////////////////
		using namespace core::scoring;
		using protocols::simple_moves::MinMoverOP;
		using protocols::simple_moves::MinMover;
		using protocols::simple_moves::TaskAwareMinMoverOP;
		using protocols::simple_moves::TaskAwareMinMover;

		ScoreFunctionOP sf(get_score_function());

		//TA + minmover for all positions fixed
		protocols::simple_moves::MinMoverOP min_mover_all( new protocols::simple_moves::MinMover(
																						mm_start,
																						sf,
																						"dfpmin_armijo", //faster and irrelevant for the purpose...
																						0.01,
																						true /*use_nblist*/ ) );
		protocols::simple_moves::TaskAwareMinMoverOP TAmin_mover_all( new protocols::simple_moves::TaskAwareMinMover(min_mover_all, task_factory_all) );

		//TA + minmover for one position fixed
		protocols::simple_moves::MinMoverOP min_mover_one( new protocols::simple_moves::MinMover(
																						mm_start,
																						sf,
																						"dfpmin_armijo", //faster and irrelevant for the purpose...
																						0.01,
																						true /*use_nblist*/ ) );
		protocols::simple_moves::TaskAwareMinMoverOP TAmin_mover_one( new protocols::simple_moves::TaskAwareMinMover(min_mover_one, task_factory_one) );

		//TA + minmover for no positions fixed
		protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover(
																				mm_start,
																				sf,
																				"dfpmin_armijo", //faster and irrelevant for the purpose...
																				0.01,
																				true /*use_nblist*/ ) );
		protocols::simple_moves::TaskAwareMinMoverOP TAmin_mover( new protocols::simple_moves::TaskAwareMinMover(min_mover, task_factory) );


		//now the not-task-aware part
		core::kinematics::MoveMapOP mm_all( new core::kinematics::MoveMap() );
		mm_all->set_bb(false); mm_all->set_jump(false); mm_all->set_chi(false);

		core::kinematics::MoveMapOP mm_one( new core::kinematics::MoveMap() );
		mm_one->set_bb(false); mm_one->set_jump(false); mm_one->set_chi(1, false); mm_one->set_chi(2, true);

		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap() );
		mm->set_bb(false); mm->set_jump(false); mm->set_chi(true);

		//MinMover for all positions fixed
		protocols::simple_moves::MinMoverOP noTA_min_mover_all( new protocols::simple_moves::MinMover(
																								 mm_all,
																								 sf,
																								 "dfpmin_armijo", //faster and irrelevant for the purpose...
																								 0.01,
																								 true /*use_nblist*/ ) );

		//MinMover for one position fixed
		protocols::simple_moves::MinMoverOP noTA_min_mover_one( new protocols::simple_moves::MinMover(
																								 mm_one,
																								 sf,
																								 "dfpmin_armijo", //faster and irrelevant for the purpose...
																								 0.01,
																								 true /*use_nblist*/ ) );

		//MinMover for no positions fixed
		protocols::simple_moves::MinMoverOP noTA_min_mover( new protocols::simple_moves::MinMover(
																						 mm,
																						 sf,
																						 "dfpmin_armijo", //faster and irrelevant for the purpose...
																						 0.01,
																						 true /*use_nblist*/ ) );

		//now make comparisons via pose coordinates compare_atom_coordinates
		core::Size const precision(6);

		//check that we're ok to start
		TS_ASSERT(core::pose::compare_atom_coordinates(p1, p2, precision));

		//compare all positions fixed
		noTA_min_mover_all->apply(p1);
		TAmin_mover_all->apply(p2);
		TS_ASSERT(core::pose::compare_atom_coordinates(p1, p2, precision));
		p1 = backup; p2 = backup;

		//compare one position fixed
		noTA_min_mover_one->apply(p1);
		TAmin_mover_one->apply(p2);
		TS_ASSERT(core::pose::compare_atom_coordinates(p1, p2, precision));
		p1 = backup; p2 = backup;

		//compare no positions fixed
		noTA_min_mover->apply(p1);
		TAmin_mover->apply(p2);
		TS_ASSERT(core::pose::compare_atom_coordinates(p1, p2, precision));
		p1 = backup; p2 = backup;

	}//end test_TaskAwareMinMover

};//end class
