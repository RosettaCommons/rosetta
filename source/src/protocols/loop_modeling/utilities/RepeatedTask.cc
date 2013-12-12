// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/utilities/RepeatedTask.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace loop_modeling {
namespace utilities {

using namespace std;

using core::Size;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using protocols::loops::Loop;

RepeatedTask::RepeatedTask(LoopMoverTaskOP task, Size iterations) {
	task_ = task;
	iterations_ = iterations;
}

void RepeatedTask::setup(
		Pose & pose, Loop const & loop, ScoreFunctionOP score_function) {

	task_->setup(pose, loop, score_function);
}

bool RepeatedTask::apply(
		Pose & pose, Loop const & loop, ScoreFunctionCOP score_function) {

	for (Size i = 1; i <= iterations_; i++) {
		bool success = task_->apply(pose, loop, score_function);
		if (success == false) return false;
	}
	return true;
}

void RepeatedTask::debrief(bool was_accepted) {
	task_->debrief(was_accepted);
}

} // namespace utilities
} // namespace kinematic_closure
} // namespace protocols

