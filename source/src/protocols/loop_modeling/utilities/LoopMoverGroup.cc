// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/utilities/LoopMoverGroup.hh>
#include <protocols/loop_modeling/utilities/LoopFilter.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/filters/Filter.hh>

// Utility headers
#include <boost/foreach.hpp>
#include <utility/vector1.hh>

#define foreach BOOST_FOREACH

namespace protocols {
namespace loop_modeling {
namespace utilities {

using namespace std;

LoopMoverGroup::LoopMoverGroup() // {{{1
: is_default_(false) {}

bool LoopMoverGroup::do_apply(Pose & pose) { // {{{1
	foreach ( LoopMoverOP child, get_children() ) {
		child->apply(pose);
		if ( ! child->was_successful() ) return false;
	}
	return true;
}

void LoopMoverGroup::get_children_names( // {{{1
	utility::vector1<string> & names, string indent) const {

	foreach ( LoopMoverOP mover, get_children() ) {
		mover->get_children_names(names, indent);
	}
}
LoopMoverOP LoopMoverGroup::add_mover(LoopMoverOP mover) { // {{{1
	if ( is_default_ ) { clear(); }
	return add_child(mover);
}

LoopMoverOP LoopMoverGroup::add_filter(protocols::filters::FilterOP filter) { // {{{1
	return add_mover(LoopMoverOP( new utilities::LoopFilter(filter) ));
}

LoopMoverGroupOP LoopMoverGroup::add_mover_group() { // {{{1
	LoopMoverGroupOP group( new LoopMoverGroup );
	add_mover(group);
	return group;
}

void LoopMoverGroup::clear() { // {{{1
	clear_children();
	is_default_ = false;
}

void LoopMoverGroup::mark_as_default() { // {{{1
	is_default_ = true;
}

bool LoopMoverGroup::empty() const { // {{{1
	return count_children() == 0;
}
// }}}1

}
}
}
