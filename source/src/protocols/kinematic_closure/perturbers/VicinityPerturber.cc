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
#include <protocols/kinematic_closure/perturbers/VicinityPerturber.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>

// Utility headers
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/kinematic_closure/vector.hh>

// C++ headers
#include <algorithm>


namespace protocols {
namespace kinematic_closure {
namespace perturbers {

using numeric::random::uniform;
using numeric::random::gaussian;
using numeric::conversions::DEGREES;

/// @details Note that the given target must have the same number of residues
/// as the pose being sampled.

VicinityPerturber::VicinityPerturber(Pose const & target)
: target_(target), spread_(10 /*degrees*/) {}

void VicinityPerturber::perturb_subset(Pose const &, IndexList const &, ClosureProblemOP problem)
{
	for ( Size const i : problem->nonpivot_residues() ) {
		problem->perturb_phi(i, target_.phi(i) + spread_ * gaussian(), DEGREES);
		problem->perturb_psi(i, target_.psi(i) + spread_ * gaussian(), DEGREES);
	}
}

void VicinityPerturber::perturb_subset_with_balance(Pose const &, IndexList const &, ClosureProblemOP problem)
{
	Real fwhm = 2. * sqrt(2. * log(2.)) * spread_;

	for ( Size const i : problem->nonpivot_residues() ) {
		Real phi = target_.phi(i) + fwhm * (uniform() - 0.5);
		Real psi = target_.psi(i) + fwhm * (uniform() - 0.5);

		problem->perturb_phi(i, phi, DEGREES);
		problem->perturb_psi(i, psi, DEGREES);
	}
}

}
}
}
