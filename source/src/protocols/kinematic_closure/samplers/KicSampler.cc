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
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>
#include <protocols/kinematic_closure/samplers/KicSampler.hh>
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
#include <numeric/random/random.hh>
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

namespace protocols {
namespace kinematic_closure {
namespace samplers {

using namespace std;
using protocols::kinematic_closure::pivot_pickers::PivotPickerOP;
using protocols::kinematic_closure::solution_pickers::SolutionPickerOP;
using protocols::loop_modeling::loggers::LoggerOP;

KicSampler::KicSampler() { // {{{1
	perturbers_ = new perturbers::PerturberSet;
	perturbers_->add(new perturbers::RamaPerturber);
	perturbers_->add(new perturbers::BondAnglePerturber);
	perturbers_->mark_as_default();

	pivot_picker_ = new pivot_pickers::StandardPivots;
	solution_picker_ = new solution_pickers::FilteredSolutions;
	logger_ = new protocols::loop_modeling::loggers::NullLogger;
	setup_called_ = false;
}

KicSampler::~KicSampler() {} // {{{1

// {{{1
/// @details This function needs to be called whenever the fold tree needs to 
/// be reset.  It is called automatically before apply() is called for the 
/// first time on this object.  After that, the user is responsible for 
/// manually calling this method.  This would be necessary, for example, if 
/// some residues were either inserted into or deleted from the structure.

void KicSampler::setup(Pose & pose, Loop const & loop) {
	setup_fold_tree(pose, loop);
	setup_called_ = true;
}

bool KicSampler::apply(Pose & pose, Loop const & loop) { // {{{1
	ClosureProblemOP problem = new ClosureProblem(logger_);
	SolutionList solutions;

	bool problem_solved = false;
	if (not setup_called_) setup(pose, loop);

	problem->frame(pose, loop, pivot_picker_);

	for (Size i = 1; i <= 2000 and not problem_solved; i++) {
		perturbers_->perturb(pose, problem);
		problem->solve(solutions);

		logger_->log_task(pose, "BridgeObjects", solutions.size() > 0, false);
		problem_solved = solution_picker_->pick_and_apply(pose, solutions);
	}

	return problem_solved;
}

// {{{1
/// @details The KicSampler starts off with a default set of perturbers.  Then 
/// first time this method is called, all of the default perturbers are cleared 
/// and only the specified perturber is kept.  Subsequent calls to this method 
/// simply add the specified perturbers.  Nothing you add manually gets erased.

void KicSampler::add_perturber(perturbers::PerturberOP perturber) {
	perturbers_->add(perturber);
}

void KicSampler::set_pivot_picker(PivotPickerOP picker) { // {{{1
	pivot_picker_ = picker;
}

void KicSampler::set_solution_picker(SolutionPickerOP picker) { // {{{1
	solution_picker_ = picker;
}

void KicSampler::log_filters(LoggerOP logger) { // {{{1
	logger_ = logger;
}
// }}}1

}
}
}
