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
#include <protocols/kinematic_closure/pivot_pickers/EndToEndPivots.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>

// Utility headers
#include <utility/exit.hh>
#include <numeric/random/random.hh>

namespace protocols {
namespace kinematic_closure {
namespace pivot_pickers {

Loop EndToEndPivots::pick(Pose const & pose, Loop const & loop) {
	runtime_assert(loop.start() < loop.stop());

	Size pivot_1 = loop.start();
	Size pivot_3 = loop.stop();
	Size pivot_2 = numeric::random::random_range(pivot_1 + 1, pivot_3 - 1);

	return Loop(pivot_1, pivot_3, pivot_2);
}

}
}
}
