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
#include <protocols/loop_modeling/utilities/LoopMoverGroup.hh>
#include <protocols/loop_modeling/utilities/LoopFilter.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/filters/Filter.hh>

// Utility headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

namespace protocols {
namespace loop_modeling {
namespace utilities {

using namespace std;
using core::scoring::ScoreFunctionOP;

bool LoopMoverGroup::do_apply(Pose & pose) { // {{{1
	foreach (LoopMoverOP child, children_) {
		child->apply(pose);
		if (not child->was_successful()) return false;
	}
	return true;
}

void LoopMoverGroup::add_mover(LoopMoverOP mover) { // {{{1
	register_nested_loop_mover(mover);
	children_.push_back(mover);
}

void LoopMoverGroup::add_filter(protocols::filters::FilterOP filter) { // {{{1
	add_mover(new utilities::LoopFilter(filter));
}

bool LoopMoverGroup::empty() const { // {{{1
	return children_.empty();
}
// }}}1

}
}
}
