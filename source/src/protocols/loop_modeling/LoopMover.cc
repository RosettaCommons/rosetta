// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/LoopMoverTask.hh>
#include <protocols/loop_modeling/utilities/FilterTask.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/filters/Filter.hh>

// Utility headers
#include <boost/foreach.hpp>

namespace protocols {
namespace loop_modeling {

#define foreach BOOST_FOREACH
using namespace std;

LoopMover::LoopMover() {
	score_function_ = NULL;
	setup_needed_ = true;
}

LoopMover::LoopMover(Loop const & loop, ScoreFunctionOP score_function) {
	set_loop(loop);
	set_score_function(score_function);
}

void LoopMover::apply(Pose & pose) {
	if (setup_needed_) setup(pose);

	last_move_successful_ = true;

	foreach (LoopMoverTaskOP task, tasks_) {
		last_move_successful_ = task->apply(pose, loop_, score_function_);

		foreach (loggers::LoggerOP logger, loggers_) {
			logger->log_task(pose, task->get_name(), last_move_successful_);
		}

		if (not last_move_successful_) break;
	}
}

void LoopMover::setup(Pose & pose) {
	if (loop_.length() == 0) {
		utility_exit_with_message("No loop region specified");
	}
	if (not score_function_) {
		score_function_ = core::scoring::getScoreFunction();
	}

	foreach (LoopMoverTaskOP task, tasks_) {
		task->setup(pose, loop_, score_function_);
	}

	setup_needed_ = false;
}

void LoopMover::set_score_function(ScoreFunctionOP score_function) {
	score_function_ = score_function;
	setup_needed_ = true;
}

void LoopMover::set_loop(Loop const & loop) {
	loop_ = loop;
	setup_needed_ = true;
}

void LoopMover::add_task(LoopMoverTaskOP task) {
	tasks_.push_back(task);
}

void LoopMover::add_filter(protocols::filters::FilterOP filter) {
	tasks_.push_back(new utilities::FilterTask(filter));
}

void LoopMover::add_logger(loggers::LoggerOP logger) {
	loggers_.push_back(logger);
}

Loop LoopMover::get_loop() const {
	return loop_;
}

ScoreFunctionOP LoopMover::get_score_function() const {
	return score_function_;
}

bool LoopMover::was_successful() const {
	return last_move_successful_;
}

}
}

