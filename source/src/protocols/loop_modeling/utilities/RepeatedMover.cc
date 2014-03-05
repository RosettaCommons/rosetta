// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/utilities/RepeatedMover.hh>

// Core headers
#include <core/pose/Pose.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

RepeatedMover::RepeatedMover(LoopMoverOP mover, Size iterations)
	: mover_(mover), iterations_(iterations) {

	register_nested_loop_mover(mover);
}

bool RepeatedMover::do_apply(Pose & pose) {
	for (Size i = 1; i <= iterations_; i++) {
		mover_->apply(pose);
		if (not mover_->was_successful()) return false;
	}
	return true;
}

} // namespace utilities
} // namespace kinematic_closure
} // namespace protocols

