// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/MinMover.cxxtest.hh
/// @brief  test for MinMover
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Project Headers
#include <core/types.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility/tag/Tag.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

// Utility Headers
#include <basic/Tracer.hh>

// C++ Headers
#include <sstream>

static basic::Tracer TR("protocols.simple_moves.MinMover.cxxtest.hh");

// --------------- Test Class --------------- //

class MyMinMover : public protocols::simple_moves::MinMover {

public:
	void
	public_apply_dof_tasks_to_movemap(
		core::pose::Pose const & pose,
		core::kinematics::MoveMap & movemap) const {
		apply_dof_tasks_to_movemap(pose, movemap);
	}
};

class MinMoverTests : public CxxTest::TestSuite {

private:
	core::pose::Pose pose_;
	core::scoring::ScoreFunctionOP scorefxn_;
public:

	void setUp() {
    core_init_with_additional_options( "-mute core.init core.pack.task core.conformation" );
    core::import_pose::pose_from_pdb( pose_, "core/pack/task/resfile_test.pdb" );
		scorefxn_ = new core::scoring::ScoreFunction;
	}

	void tearDown() {
	}

	void test_use_setup_dofs_from_resfile() {
		using namespace protocols::simple_moves;
		using core::pack::task::operation::TaskOperationOP;
		using core::pack::task::operation::ReadResfile;
		using namespace utility::tag;
		using core::scoring::ScoreFunctionFactory;

		core::kinematics::MoveMapOP ref_movemap(new core::kinematics::MoveMap());
		ref_movemap->set_bb(true);
		ref_movemap->set_bb(3, false);
		ref_movemap->set_bb(6, false);
		ref_movemap->set_bb(9, false);
		ref_movemap->set_bb(10, false);
		ref_movemap->set_bb(12, false);
		ref_movemap->set_bb(13, false);
		ref_movemap->set_bb(14, false);
		ref_movemap->set_bb(16, false);
		ref_movemap->set_bb(17, false);
		ref_movemap->set_bb(18, false);
		ref_movemap->set_bb(19, false);
		ref_movemap->set_bb(20, false);
		ref_movemap->set_bb(27, false);
		ref_movemap->set_bb(28, false);
		ref_movemap->set_bb(29, false);
		ref_movemap->set_bb(30, false);
		TR << "Reference MoveMap" << std::endl;
		TR << *ref_movemap << std::endl;

		MyMinMover min_mover;

		protocols::moves::DataMap data_map;

		data_map.add(
			"task_operations", "resfile",
			TaskOperationOP( new ReadResfile("core/pack/task/resfile_test.resfile")));

		data_map.add(
			"scorefxns", "score12",
			core::scoring::ScoreFunctionFactory::create_score_function("score12"));

		std::stringstream tag_st("<MinMover chi=0 bb=1 bb_task_operations=resfile/>");
		TagPtr tag = Tag::create(tag_st);

		protocols::filters::Filters_map filters;
		protocols::moves::Movers_map movers;

		min_mover.parse_my_tag(tag, data_map, filters, movers, pose_);
		core::kinematics::MoveMapOP movemap(min_mover.movemap()->clone());
		min_mover.public_apply_dof_tasks_to_movemap(pose_, *movemap);

		TR << "After apply_dof_tasks_to_movemap" << std::endl;
		TR << *movemap << std::endl;

		TS_ASSERT(*ref_movemap == *movemap);
	}

};
