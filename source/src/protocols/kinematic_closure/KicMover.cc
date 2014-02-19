// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/utilities.hh>
#include <protocols/kinematic_closure/KicMover.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/perturbers/PerturberSet.hh>
#include <protocols/kinematic_closure/perturbers/RamaPerturber.hh>
#include <protocols/kinematic_closure/perturbers/BondAnglePerturber.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.hh>
#include <protocols/kinematic_closure/pivot_pickers/StandardPivots.hh>
#include <protocols/kinematic_closure/solution_pickers/SolutionPicker.hh>
#include <protocols/kinematic_closure/solution_pickers/FilteredSolutions.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/loop_modeling/loggers/Logger.hh>
#include <protocols/loop_modeling/loggers/NullLogger.hh>

// Utility headers
#include <utility/exit.hh>
#include <numeric/random/random.hh>
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

namespace protocols {
namespace kinematic_closure {

using namespace std;
using protocols::kinematic_closure::pivot_pickers::PivotPickerOP;
using protocols::kinematic_closure::solution_pickers::SolutionPickerOP;
using protocols::loop_modeling::loggers::LoggerOP;

KicMover::KicMover() { // {{{1
	perturbers_ = new perturbers::PerturberSet;
	perturbers_->add(new perturbers::RamaPerturber);
	perturbers_->add(new perturbers::BondAnglePerturber);
	perturbers_->mark_as_default();

	pivot_picker_ = new pivot_pickers::StandardPivots;
	solution_picker_ = new solution_pickers::FilteredSolutions;
	logger_ = new protocols::loop_modeling::loggers::NullLogger;
	is_fold_tree_stale_ = true;
}

KicMover::~KicMover() {} // {{{1

void KicMover::apply(Pose & pose) { // {{{1
	ClosureProblemOP problem = new ClosureProblem(logger_);
	SolutionList solutions;

	if (loop_.length() == 0) {
		utility_exit_with_message(
				"Before calling BalancedKicMover.apply(), you must provide a loop "
				"via BalancedKicMover.set_loop().");
	}

	bool problem_solved = false;
	problem->frame(pose, loop_, pivot_picker_);

	// Make 2000 attempts to find a closure solution which passes which passes 
	// both rama and bump checks.  If no solution is found by then, give up.

	for (Size i = 1; i <= 2000 and not problem_solved; i++) {
		perturbers_->perturb(pose, problem);
		problem->solve(solutions);

		logger_->log_task(pose, "BridgeObjects", solutions.size() > 0, false);
		problem_solved = solution_picker_->pick_and_apply(pose, solutions);
	}

	// Set the move type based on whether or not the move succeeded.
	
	type(problem_solved ? "kic" : "kic (no-op)");
}

void KicMover::setup(Pose & pose, Loop const & loop) { // {{{1
	if (loop_ != loop) {
		loop_ = loop;
		setup_fold_tree(pose, loop_);
	}
}

// {{{1
/// @details The KicMover starts off with a default set of perturbers.  Then 
/// first time this method is called, all of the default perturbers are cleared 
/// and only the specified perturber is kept.  Subsequent calls to this method 
/// simply add the specified perturbers.  Nothing you add manually gets erased.

void KicMover::add_perturber(perturbers::PerturberOP perturber) {
	perturbers_->add(perturber);
}

void KicMover::set_pivot_picker(PivotPickerOP picker) { // {{{1
	pivot_picker_ = picker;
}

void KicMover::set_solution_picker(SolutionPickerOP picker) { // {{{1
	solution_picker_ = picker;
}

void KicMover::log_filters(LoggerOP logger) { // {{{1
	logger_ = logger;
}
// }}}1

}
}
