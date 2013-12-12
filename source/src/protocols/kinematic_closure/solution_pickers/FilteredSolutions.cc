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
#include <protocols/kinematic_closure/solution_pickers/FilteredSolutions.hh>

// Core headers
#include <core/pose/Pose.hh>

// Utility headers
#include <numeric/random/random.hh>
#include <boost/foreach.hpp>

// C++ headers
#include <iostream>
#include <sstream>

#define foreach BOOST_FOREACH

namespace protocols {
namespace kinematic_closure {
namespace solution_pickers {

using namespace std;

FilteredSolutions::FilteredSolutions(bool check_rama, bool check_overlap) {
	check_rama_ = check_rama;
	check_overlap_ = check_overlap;
}

bool FilteredSolutions::pick_and_apply(
		Pose & pose, SolutionList const & solutions) {

	foreach (ClosureSolutionCOP solution, solutions) {
		bool reasonable_solution =
			solution->apply_if_reasonable(pose, check_rama_, check_overlap_);

		if (reasonable_solution) return true;
	}

	return false;
}

}
}
}
