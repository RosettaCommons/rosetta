// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/pivot_pickers/FixedOffsetPivots.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>

// Utility headers
#include <utility/exit.hh>
#include <numeric/random/random.hh>

// C++ headers
#include <algorithm>
#include <iostream>

namespace protocols {
namespace kinematic_closure {
namespace pivot_pickers {

using namespace std;

Loop FixedOffsetPivots::pick(Pose const &, Loop const & loop) {
	using numeric::random::random_range;

	Size pivot_1 = random_range(loop.start(), loop.stop() - offset_);
	Size pivot_3 = pivot_1 + offset_;
	Size pivot_2 = pivot_1 + offset_ / 2;

	return Loop(pivot_1, pivot_3, pivot_2);
}

}
}
}
