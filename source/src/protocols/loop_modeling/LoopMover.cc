// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Headers {{{1
#include <protocols/loop_modeling/LoopMover.hh>

// Core headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/util/kinematics_util.hh>

// Protocol headers
#include <protocols/loops/loops_main.hh>
#include <protocols/loop_modeling/utilities/rosetta_scripts.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/HierarchicalDataMap.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// RosettaScripts headers
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>

// }}}1

namespace protocols {
namespace loop_modeling {

// Global Names {{{1
using namespace std;
using core::kinematics::FoldTree;

static THREAD_LOCAL basic::Tracer TR( "protocols.loop_modeling.LoopMover" );
const string ToolboxKeys::LOOPS = "loops";
const string ToolboxKeys::SCOREFXN = "scorefxn";
const string ToolboxKeys::TASK_FACTORY = "task_factory";
// }}}1

LoopMover::LoopMover()  // {{{1
: children_(),
	toolbox_(new basic::datacache::HierarchicalDataMap),
	parent_name_(""),
	trust_fold_tree_(false),
	was_successful_(false) {}

void LoopMover::parse_my_tag( // {{{1
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &) {

	//SML Jul 11 2016, moved this function to a utilities section for sharing
	LoopsOP parsed_loops(protocols::loop_modeling::utilities::parse_loops_from_tag(tag));

	// Don't override any loops that may be specified in parent movers unless
	// loops really were specified here.

	if ( parsed_loops ) {
		set_loops(parsed_loops);
	}
}

void LoopMover::apply(Pose & pose) { // {{{1

	// If the existing fold tree is trusted, then just apply the loop mover.
	// Otherwise, construct an compatible fold tree, apply the loop mover, then
	// restore the existing fold tree.

	if ( trust_fold_tree_ ) {
		was_successful_ = do_apply(pose);
	} else {
		FoldTree original_tree = pose.fold_tree();
		FoldTreeRequest request = request_fold_tree();
		setup_fold_tree(pose, get_loops(), request);
		protocols::loops::add_cutpoint_variants(pose);

		was_successful_ = do_apply(pose);

		// There is a protocols::loops::remove_cutpoint_variants() function, but it
		// clobbers the constraint weights in the pose (even though it seems
		// explicitly designed to not do that).  Because this leads to incorrect
		// scores being reported in the output pdb file, I chose to use the
		// core::util function instead, which doesn't have this weakness.

		core::util::remove_cutpoint_variants(pose);
		pose.fold_tree(original_tree);
	}
}

bool LoopMover::do_apply(Pose & pose) { // {{{1
	LoopsCOP loops = get_loops();
	Loops::const_iterator loop, end;

	for ( loop = loops->begin(), end = loops->end(); loop != end; ++loop ) {
		bool was_successful = do_apply(pose, *loop);
		if ( ! was_successful ) return false;
	}
	return true;
}

bool LoopMover::do_apply(Pose &, Loop const &) { // {{{1
	utility_exit_with_message("LoopMover::do_apply was not reimplemented.");
	return false;
}

void LoopMover::get_children_names( // {{{1
	utility::vector1<string> & names, string indent) const {

	names.push_back(indent + get_name());
	for ( LoopMoverOP mover : get_children() ) {
		mover->get_children_names(names, indent + "  ");
	}
}

bool LoopMover::was_successful() const { // {{{1
	return was_successful_;
}

LoopsOP LoopMover::get_loops() { // {{{1
	return get_tool<LoopsOP>(ToolboxKeys::LOOPS);
}

LoopsCOP LoopMover::get_loops() const { // {{{1
	return get_tool<LoopsCOP>(ToolboxKeys::LOOPS);
}

Loop const & LoopMover::get_loop(Size index) const { // {{{1
	return (*get_loops())[index];
}

void LoopMover::set_loops(LoopsOP loops) { // {{{1
	set_tool(ToolboxKeys::LOOPS, loops);
}

void LoopMover::set_loops(Loops const & loops) { // {{{1
	LoopsOP loops_op( new Loops(loops) );
	set_loops(loops_op);
}

void LoopMover::set_loop(Loop const & loop) { // {{{1
	LoopsOP wrapper( new Loops );
	wrapper->add_loop(loop);
	set_loops(wrapper);
}

FoldTreeRequest LoopMover::request_fold_tree() const { // {{{1
	FoldTreeRequest request = FTR_DONT_CARE;
	for ( LoopMoverCOP mover : get_children() ) {
		request = request & mover->request_fold_tree();
	}
	return request;
}

void LoopMover::trust_fold_tree() { // {{{1
	trust_fold_tree_ = true;
}

void LoopMover::setup_fold_tree( // {{{1
	Pose & pose, LoopsCOP loops, FoldTreeRequest request) {

	if ( request == FTR_DONT_CARE ) {
		return;
	} else if ( request & FTR_LOOPS_WITH_CUTS ) {
		core::kinematics::FoldTree tree;
		protocols::loops::fold_tree_from_loops(pose, *loops, tree, true);
		pose.fold_tree(tree);
	} else if ( request & FTR_SIMPLE_TREE ) {
		core::kinematics::FoldTree tree(pose.size());
		pose.fold_tree(tree);
	} else {
		utility_exit_with_message("Could not setup the requested fold tree.");
	}

	protocols::loops::add_cutpoint_variants(pose);
}

void LoopMover::remove_child(LoopMoverOP child) { // {{{1
	child->parent_name_ = "";
	child->toolbox_->unset_parent();
	children_.erase(std::find(children_.begin(), children_.end(), child));
}

void LoopMover::clear_children() { // {{{1
	for ( LoopMoverOP child : children_ ) {
		child->parent_name_ = "";
		child->toolbox_->unset_parent();
	}
	children_.clear();
}

LoopMoverOPs LoopMover::get_children() const { // {{{1
	return children_;
}

Size LoopMover::count_children() const { // {{{1
	return children_.size();
}
// }}}1

}
}
