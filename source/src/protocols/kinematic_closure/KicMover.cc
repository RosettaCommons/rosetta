// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Headers {{{1
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/internal.hh>
#include <protocols/kinematic_closure/KicMover.hh>
#include <protocols/kinematic_closure/KicMoverCreator.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/ClosureSolution.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>
#include <protocols/kinematic_closure/perturbers/PerturberSet.hh>
#include <protocols/kinematic_closure/perturbers/Rama2bPerturber.hh>
#include <protocols/kinematic_closure/perturbers/OmegaPerturber.hh>
#include <protocols/kinematic_closure/pivot_pickers/PivotPicker.hh>
#include <protocols/kinematic_closure/pivot_pickers/StandardPivots.hh>
#include <protocols/kinematic_closure/solution_pickers/SolutionPicker.hh>
#include <protocols/kinematic_closure/solution_pickers/FilteredSolutions.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Basic headers
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/option.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/moves/Mover.hh>

// Utility headers
#include <boost/foreach.hpp>
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// Namespaces {{{1
using namespace std;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using protocols::kinematic_closure::pivot_pickers::PivotPickerOP;
using protocols::kinematic_closure::solution_pickers::SolutionPickerOP;
using protocols::loop_modeling::FoldTreeRequest;
using protocols::loop_modeling::FTR_LOOPS_WITH_CUTS;
// }}}1

namespace protocols {
namespace kinematic_closure {

protocols::moves::MoverOP KicMoverCreator::create_mover() const { // {{{1
	return protocols::moves::MoverOP( new KicMover );
}

std::string KicMoverCreator::keyname() const { // {{{1
	return "KicMover";
}
// }}}1

KicMover::KicMover() { // {{{1
	using perturbers::PerturberOP;

	perturbers_ = perturbers::PerturberSetOP( new perturbers::PerturberSet );
	perturbers_->add(PerturberOP( new perturbers::Rama2bPerturber ));
	perturbers_->add(PerturberOP( new perturbers::OmegaPerturber ));
	perturbers_->mark_as_default();

	pivot_picker_ = pivot_pickers::PivotPickerOP( new pivot_pickers::StandardPivots );
	solution_picker_ = solution_pickers::SolutionPickerOP( new solution_pickers::FilteredSolutions );
}

KicMover::~KicMover() {} // {{{1

bool KicMover::do_apply(Pose & pose, Loop const & loop) { // {{{1
	ClosureProblemOP problem( new ClosureProblem() );
	problem->frame(pose, loop, pivot_picker_);

	bool problem_solved = false;
	Size const max_attempts = option[OptionKeys::loops::max_kic_perturber_samples]();

	// Attempt to find a closure solution which passes both rama and bump checks.
	// If no solution is found, give up.

	for ( Size i = 1; i <= max_attempts && !problem_solved; i++ ) {
		perturbers_->perturb(pose, problem);
		SolutionList solutions = problem->solve();
		problem_solved = solution_picker_->pick_and_apply(pose, solutions);
	}

	// Set the move type based on whether or not the move succeeded.

	type(problem_solved ? "kic" : "kic-no-op");
	return problem_solved;
}

void KicMover::get_children_names( // {{{1
	utility::vector1<string> & names, string indent) const {

	names.push_back(indent + get_name());
	names.push_back(indent + "  " + pivot_picker_->get_name());
	perturbers_->get_perturber_names(names, indent + "  ");
}

pivot_pickers::PivotPickerOP KicMover::get_pivot_picker() { // {{{1
	return pivot_picker_;
}

solution_pickers::SolutionPickerOP KicMover::get_solution_picker() { // {{{1
	return solution_picker_;
}

// {{{1
/// @details The KicMover starts off with a default set of perturbers.  Then
/// first time this method is called, all of the default perturbers are cleared
/// and only the specified perturber is kept.  Subsequent calls to this method
/// simply add the specified perturbers.  Nothing you add manually gets erased.

void KicMover::add_perturber(perturbers::PerturberOP perturber) {
	perturbers_->add(perturber);
}

void KicMover::clear_perturbers() { // {{{1
	perturbers_->clear();
}

void KicMover::set_pivot_picker(PivotPickerOP picker) { // {{{1
	pivot_picker_ = picker;
}

void KicMover::set_solution_picker(SolutionPickerOP picker) { // {{{1
	solution_picker_ = picker;
}

FoldTreeRequest KicMover::request_fold_tree() const { // {{{1
	return FTR_LOOPS_WITH_CUTS;
}
// }}}1

}
}
