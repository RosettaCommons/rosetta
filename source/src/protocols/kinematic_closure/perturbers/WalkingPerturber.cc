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
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/perturbers/WalkingPerturber.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>

// Utility headers
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/kinematic_closure/vector.hh>


namespace protocols {
namespace kinematic_closure {
namespace perturbers {

using numeric::random::gaussian;
using numeric::conversions::DEGREES;

WalkingPerturber::WalkingPerturber(Real magnitude)
: magnitude_ (magnitude) {}

void WalkingPerturber::perturb_subset(
	Pose const &, IndexList const & residues, ClosureProblemOP problem) {

	for ( Size const residue : residues ) {
		Real phi = problem->phi(residue, DEGREES) + magnitude_ * gaussian();
		Real psi = problem->psi(residue, DEGREES) + magnitude_ * gaussian();

		problem->perturb_phi(residue, phi, DEGREES);
		problem->perturb_psi(residue, psi, DEGREES);
	}
}

void WalkingPerturber::perturb_subset_with_balance(
	Pose const & pose, IndexList const & residues, ClosureProblemOP problem) {

	// Each move makes a gaussian step in a random direction.  Since the choice
	// of direction is unbiased, the overall move satisfies detailed balance.

	perturb_subset(pose, residues, problem);
}

}
}
}
