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
#include <protocols/kinematic_closure/perturbers/WalkingBondAnglePerturber.hh>

// Core headers
#include <core/pose/Pose.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/kinematic_closure/vector.hh>
#include <boost/foreach.hpp>


namespace protocols {
namespace kinematic_closure {
namespace perturbers {

WalkingBondAnglePerturber::WalkingBondAnglePerturber(Real magnitude)
: magnitude_ (magnitude) {}

void WalkingBondAnglePerturber::perturb_subset(Pose const &, IndexList const & residues, ClosureProblemOP problem) {
	using numeric::random::gaussian;
	using numeric::conversions::DEGREES;

	BOOST_FOREACH ( Size residue, residues ) {
		Real angle = problem->n_ca_c(residue, DEGREES) + magnitude_ * gaussian();
		problem->perturb_n_ca_c(residue, angle, DEGREES);
	}
}

void WalkingBondAnglePerturber::perturb_subset_with_balance(
	Pose const & pose, IndexList const & residues, ClosureProblemOP problem) {

	// Each move makes a gaussian step in a random direction.  Since the choice
	// of direction is unbiased, the overall move satisfies detailed balance.

	perturb_subset(pose, residues, problem);
}

}
}
}
