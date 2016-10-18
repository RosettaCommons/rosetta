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
#include <protocols/kinematic_closure/ClosureSolution.hh>
#include <protocols/kinematic_closure/solution_pickers/FilteredSolutions.hh>

// Core headers
#include <core/pose/Pose.hh>

// Utility headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <boost/foreach.hpp>

// C++ headers
#include <iostream>
#include <sstream>


namespace protocols {
namespace kinematic_closure {
namespace solution_pickers {

using namespace std;

FilteredSolutions::FilteredSolutions(
	bool check_rama, bool check_overlap, bool be_lenient) {

	check_rama_ = check_rama;
	check_overlap_ = check_overlap;
	be_lenient_ = be_lenient;
}

bool FilteredSolutions::pick_and_apply(
	Pose & pose, SolutionList const & solutions) {

	core::Size nsol = solutions.size();

	// Permute the solutions randomly before choosing one,
	// since the first one encountered that passes the filters is accepted.

	utility::vector1<core::Size> pos(nsol);
	for ( core::Size i=1; i<= nsol; ++i ) {
		pos[i] = i;
	}

	numeric::random::random_permutation(pos.begin(), pos.end(), numeric::random::rg());

	for ( core::Size i=1; i <= nsol; ++i ) {
		bool reasonable_solution =
			solutions[pos[i]]->apply_if_reasonable(
			pose, check_rama_, check_overlap_, be_lenient_);

		if ( reasonable_solution ) return true;
	}

	return false;
}

}
}
}
