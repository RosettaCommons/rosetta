// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/LoopMover.hh>

// Core headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Protocol headers
#include <protocols/loops/loops_main.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <boost/foreach.hpp>
#include <utility/exit.hh>

namespace protocols {
namespace loop_modeling {

using namespace std;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using core::kinematics::FoldTree;

static thread_local basic::Tracer TR( "protocols.loop_modeling.LoopMover" );

LoopMover::LoopMover()  // {{{1
	: trust_fold_tree_(false),
		manage_score_function_(true),
		was_successful_(false),
		score_function_(NULL)
	{}

LoopMover::~LoopMover() {}  // {{{1
// }}}1

void LoopMover::apply(Pose & pose) { // {{{1

	// Make sure a loop was specified.

	if (loops_.empty()) {
		utility_exit_with_message("No loops specified.");
	}

	// Create a score function if necessary.
	if (score_function_.get() == NULL) {
		set_score_function(core::scoring::get_score_function());
	}

	// Sample the loops.  If the existing fold tree isn't trusted, replace it
	// with a custom-built one first.

	if (trust_fold_tree_) {
		was_successful_ = do_apply(pose);
	}
	else {
		FoldTree original_tree = pose.fold_tree();
		FoldTreeRequest request = request_fold_tree();
		setup_fold_tree(pose, loops_, request);
		protocols::loops::add_cutpoint_variants(pose);

		was_successful_ = do_apply(pose);

		protocols::loops::remove_cutpoint_variants(pose);
		pose.fold_tree(original_tree);
	}
}

bool LoopMover::do_apply(Pose & pose) { // {{{1
	Loops::const_iterator loop;
	for (loop = loops_.begin(); loop != loops_.end(); loop++) {
		bool was_successful = do_apply(pose, *loop);
		if (! was_successful) return false;
	}
	return true;
}

bool LoopMover::do_apply(Pose &, Loop const &) { // {{{1
	utility_exit_with_message("LoopMover::do_apply was not reimplemented.");
	return false;
}
// }}}1

bool LoopMover::was_successful() const { // {{{1
	return was_successful_;
}

Loops LoopMover::get_loops() const { // {{{1
	return loops_;
}

ScoreFunctionCOP LoopMover::get_score_function() const { // {{{1
	if (! manage_score_function_) {
		TR.Warning << "get_score_function() being called on LoopMover which doesn't support it." << endl;
		return NULL;
	}
	return score_function_;
}

ScoreFunctionOP LoopMover::get_score_function() { // {{{1
	if (! manage_score_function_) {
		TR.Warning << "get_score_function() being called on LoopMover which doesn't support it." << endl;
		return NULL;
	}
	return score_function_;
}

void LoopMover::set_loops(Loops const & loops) { // {{{1
	loops_ = loops;
	BOOST_FOREACH (LoopMoverOP mover, nested_movers_) {
		mover->set_loops(loops);
	}
}

void LoopMover::set_loop(Loop const & loop) { // {{{1
	Loops wrapper;
	wrapper.add_loop(loop);
	set_loops(wrapper);
}

void LoopMover::set_score_function(ScoreFunctionOP function) { // {{{1
	if (! manage_score_function_) {
		TR.Warning << "set_score_function() being called on LoopMover which doesn't support it." << endl;
		return;
	}

	score_function_ = function;
	BOOST_FOREACH (LoopMoverOP mover, nested_movers_) {
		mover->set_score_function(function);
	}
}

FoldTreeRequest LoopMover::request_fold_tree() const { // {{{1
	FoldTreeRequest request = FTR_DONT_CARE;
	BOOST_FOREACH (LoopMoverCOP mover, nested_movers_) {
		request = request & mover->request_fold_tree();
	}
	return request;
}

void LoopMover::trust_fold_tree() { // {{{1
	trust_fold_tree_ = true;
}

void LoopMover::setup_fold_tree( // {{{1
		Pose & pose, Loops const & loops, FoldTreeRequest request) {

	if (request == FTR_DONT_CARE) {
		return;
	}

	else if (request & FTR_LOOPS_WITH_CUTS) {
		core::kinematics::FoldTree tree;
		protocols::loops::fold_tree_from_loops(pose, loops, tree, true);
		pose.fold_tree(tree);
	}

	else if (request & FTR_SIMPLE_TREE) {
		core::kinematics::FoldTree tree(pose.total_residue());
		pose.fold_tree(tree);
	}

	else {
		utility_exit_with_message("Could not setup the requested fold tree.");
	}

	protocols::loops::add_cutpoint_variants(pose);
}

void LoopMover::dont_manage_score_function() { // {{{1
	manage_score_function_ = false;
}

void LoopMover::deregister_nested_loop_movers() { // {{{1
	nested_movers_.clear();
}

vector1<LoopMoverOP> const & LoopMover::get_nested_loop_movers() const { // {{{1
	return nested_movers_;
}
// }}}1

}
}
