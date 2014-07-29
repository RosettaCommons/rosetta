// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/utilities/PeriodicMover.hh>

// Core headers
#include <core/pose/Pose.hh>

namespace protocols {
namespace loop_modeling {
namespace utilities {

using namespace std;

using core::Size;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using protocols::loops::Loop;

PeriodicMover::PeriodicMover(LoopMoverOP mover, Size period) {
	mover_ = register_nested_loop_mover(mover);
	period_ = period;
	iteration_ = 0;
}

bool PeriodicMover::do_apply(Pose & pose) {
	if (iteration_++ % period_ == 0) {
		mover_->apply(pose);
		return mover_->was_successful();
	}
	else return true;
}

} // namespace utilities
} // namespace loop_modeling
} // namespace protocols

