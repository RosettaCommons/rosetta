// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/openers/OpenerCollection.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>

// Utility headers
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

namespace protocols {
namespace loop_modeling {
namespace openers {

OpenerCollection::OpenerCollection() {
	is_default_ = false;
}

OpenerCollection::~OpenerCollection() {}

void OpenerCollection::apply(Pose & pose, Loop const & loop) {
	foreach (LoopOpenerOP opener, openers_) {
		opener->apply(pose, loop);
	}
}

void OpenerCollection::add(LoopOpenerOP opener) {
	if (is_default_) { clear(); is_default_ = false; }
	openers_.push_back(opener);
}

void OpenerCollection::clear() {
	openers_.clear();
}

void OpenerCollection::mark_as_default() {
	is_default_ = true;
}

}
}
}


