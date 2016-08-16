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
#include <protocols/kinematic_closure/perturbers/UniformPerturber.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>

// Utility headers
#include <numeric/constants.hh>
#include <numeric/random/random.hh>
#include <numeric/kinematic_closure/vector.hh>
#include <boost/foreach.hpp>


namespace protocols {
namespace kinematic_closure {
namespace perturbers {

void UniformPerturber::perturb_subset(
	Pose const &, IndexList const & residues, ClosureProblemOP problem) {

	using numeric::random::uniform;
	using numeric::conversions::DEGREES;

	BOOST_FOREACH ( Real residue, residues ) {
		problem->perturb_phi(residue, 360 * uniform(), DEGREES);
		problem->perturb_psi(residue, 360 * uniform(), DEGREES);
		problem->perturb_omega(residue, 360 * uniform(), DEGREES);
	}
}

void UniformPerturber::perturb_subset_with_balance(
	Pose const & pose, IndexList const & residues, ClosureProblemOP problem) {

	// Uniform moves are always balanced.  To put it more technically, the
	// forward and reverse proposal probabilities are equal a priori when moves
	// are being picked from a uniform distribution.

	perturb_subset(pose, residues, problem);
}

}
}
}
