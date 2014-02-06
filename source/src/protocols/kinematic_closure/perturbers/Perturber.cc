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
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.hh>

// Core headers
#include <core/pose/Pose.hh>

namespace protocols {
namespace kinematic_closure {
namespace perturbers {

using namespace std;

void Perturber::perturb(
		Pose const & pose, ClosureProblemOP problem) {

	perturb_subset(pose, problem->residues(), problem);
}

/// @details If not reimplemented, this method will kill rosetta and complain 
/// that no balanced implementation of this algorithm exists.

void Perturber::perturb_with_balance(
		Pose const & pose, ClosureProblemOP problem) {

	perturb_subset_with_balance(pose, problem->residues(), problem);
}

/// @details If not reimplemented, this method will kill rosetta and complain 
/// that no balanced implementation of this algorithm exists.

void Perturber::perturb_subset_with_balance(
		Pose const &, IndexList const &, ClosureProblemOP) {

	utility_exit_with_message(
			"The " + get_name() + " perturber can not sample without bias.");
}

}
}
}

