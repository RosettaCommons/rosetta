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
#include <protocols/kinematic_closure/pivot_pickers/StandardPivots.hh>

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

Loop StandardPivots::pick(Pose const &, Loop const & loop) {
	int loop_start = loop.start();
	int loop_stop = loop.stop();
	int pivot_1, pivot_2, pivot_3;

	int const minimum_size = 3;
	int const maximum_size = 12;
	int const min_offset = minimum_size - 1;
	int const max_offset = maximum_size - 1;

	using numeric::random::random_range;

	// I don't think this will give a uniform distribution of pivots.  The first
	// pivot will be uniform, but the second pivot will be more likely to be one
	// of the last few residues because of the min() call.  Flipping the
	// directional bias every time will fudge things a bit, but there will still
	// be a bias away from the middle.

	if ( counter_++ % 2 == 0 ) {
		pivot_1 = random_range(loop_start, loop_stop - min_offset);
		pivot_3 = random_range(pivot_1 + min_offset, loop_stop);
		pivot_3 = min(pivot_3, pivot_1 + max_offset);
	} else {
		pivot_3 = random_range(loop_start + min_offset, loop_stop);
		pivot_1 = random_range(loop_start, pivot_3 - min_offset);
		pivot_1 = max(pivot_1, pivot_3 - max_offset);
	}

	pivot_2 = pivot_1 + (pivot_3 - pivot_1) / 2;

	runtime_assert(pivot_1 >= loop_start);
	runtime_assert(pivot_2 > pivot_1);
	runtime_assert(pivot_3 > pivot_2);
	runtime_assert(loop_stop >= pivot_3);

	return Loop(pivot_1, pivot_3, pivot_2);
}

}
}
}
