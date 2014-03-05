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
#include <protocols/kinematic_closure/pivot_pickers/LoopPivots.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>

// C++ headers
#include <iostream>

namespace protocols {
namespace kinematic_closure {
namespace pivot_pickers {

using namespace std;

Loop LoopPivots::pick(Pose const &, Loop const & loop) {
	Loop pivots = loop;
	if (pivots.cut() <= pivots.start() or pivots.cut() >= pivots.stop()) {
		pivots.set_cut(pivots.midpoint() - 1);
	}
	return pivots;
}

}
}
}
