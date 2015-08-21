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
#include <protocols/kinematic_closure/ClosureSolution.hh>
#include <protocols/kinematic_closure/solution_pickers/RandomSolutions.hh>

// Core headers
#include <core/pose/Pose.hh>

// Utility headers
#include <numeric/random/random.hh>
#include <boost/foreach.hpp>

// C++ headers
#include <iostream>
#include <sstream>


namespace protocols {
namespace kinematic_closure {
namespace solution_pickers {

using namespace std;
using numeric::random::random_range;

bool RandomSolutions::pick_and_apply(
	Pose & pose, SolutionList const & solutions) {

	if ( solutions.empty() ) return false;

	Size index = random_range(1, solutions.size());
	solutions[index]->apply(pose);
	return true;
}

}
}
}
