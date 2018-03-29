// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Kale Kundert (kale.kundert@ucsf.edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>
#include <test/protocols/task_operations/utilities.hh>

// Unit headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/simple_task_operations/RestrictToLoops.hh>

// Core headers
#include <core/types.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>

// C++ headers
#include <iostream>
#include <string>

using namespace std;
using namespace core;
using namespace core::pack::task;
using namespace protocols::loops;

using core::pose::Pose;
using protocols::loops::Loops;
using protocols::loops::LoopsOP;
using protocols::simple_task_operations::RestrictToLoops;
using protocols::simple_task_operations::RestrictToLoopsOP;

class RestrictToLoopsTests : public CxxTest::TestSuite {

public:

	void setUp() {
		protocols_init();
		core::pose::make_pose_from_sequence(pose_, "AAAAAAAAA", "fa_standard");
		loops_ = LoopsOP( new Loops ); loops_->add_loop(2, 4); loops_->add_loop(6, 8);
	}

	void test_copy_constructor() {
		RestrictToLoopsOP restrict_to_loops( new RestrictToLoops );
		restrict_to_loops->set_design_loop(true);
		restrict_to_loops->set_loops(loops_);

		RestrictToLoopsOP clone( new RestrictToLoops(*restrict_to_loops) );
		restrict_to_loops->set_design_loop(false);
		restrict_to_loops->set_loops(NULL);

		TS_ASSERT(clone->design_loop() == true);
		TS_ASSERT((*clone->loops())[1].start() == 2);
		TS_ASSERT((*clone->loops())[1].stop() == 4);
		TS_ASSERT((*clone->loops())[2].start() == 6);
		TS_ASSERT((*clone->loops())[2].stop() == 8);
	}

	void test_set_loops() {
		RestrictToLoopsOP restrict_to_loops( new RestrictToLoops );
		restrict_to_loops->set_loops(loops_);

		test_task_operation(restrict_to_loops, pose_, "011101110", "000000000");
	}

	void test_set_loops_from_file() {
		RestrictToLoopsOP restrict_to_loops( new RestrictToLoops );
		restrict_to_loops->set_loops_from_file(
			"protocols/task_operations/loops.txt");

		test_task_operation(restrict_to_loops, pose_, "011101110", "000000000");
	}

	void test_set_design_loops() {
		RestrictToLoopsOP restrict_to_loops( new RestrictToLoops );
		restrict_to_loops->set_design_loop(true);
		restrict_to_loops->set_loops(loops_);

		test_task_operation(restrict_to_loops, pose_, "011101110", "011101110");
	}

	void test_restrict_design_to_loops() {
		RestrictToLoopsOP restrict_to_loops( new RestrictToLoops );
		restrict_to_loops->set_restrict_only_design_to_loops(true);
		restrict_to_loops->set_loops(loops_);
		//std::cout << restrict_to_loops->restrict_only_design_to_loops() << std::endl;
		test_task_operation(restrict_to_loops, pose_, "111111111", "011101110");
	}
public:

	Pose pose_;
	LoopsOP loops_;

};


