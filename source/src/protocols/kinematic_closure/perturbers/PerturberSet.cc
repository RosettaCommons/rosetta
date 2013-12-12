// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/perturbers/PerturberSet.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>

// Utility headers
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

namespace protocols {
namespace kinematic_closure {
namespace perturbers {

PerturberSet::PerturberSet() {
	is_default_ = false;
}

PerturberSet::~PerturberSet() {}

void PerturberSet::perturb_subset(
		Pose const & pose, IndexList const & residues, ClosureProblemOP problem) {

	foreach (PerturberOP perturber, perturbers_) {
		perturber->perturb_subset(pose, residues, problem);
	}
}

void PerturberSet::perturb_subset_with_balance(
		Pose const & pose, IndexList const & residues, ClosureProblemOP problem) {

	foreach (PerturberOP perturber, perturbers_) {
		perturber->perturb_subset_with_balance(pose, residues, problem);
	}
}

void PerturberSet::add(PerturberOP perturber) {
	if (is_default_) { clear(); is_default_ = false; }
	perturbers_.push_back(perturber);
}

void PerturberSet::clear() {
	perturbers_.clear();
}

void PerturberSet::mark_as_default() {
	is_default_ = true;
}

}
}
}


