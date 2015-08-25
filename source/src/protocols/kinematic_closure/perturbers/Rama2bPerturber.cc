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
#include <protocols/kinematic_closure/perturbers/Rama2bPerturber.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/Ramachandran2B.hh>
#include <core/scoring/ScoringManager.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>

// Utility headers
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/kinematic_closure/vector.hh>
#include <boost/foreach.hpp>


namespace protocols {
namespace kinematic_closure {
namespace perturbers {

using core::scoring::Ramachandran2B;
using core::scoring::ScoringManager;
using numeric::conversions::DEGREES;

void Rama2bPerturber::perturb_subset(
	Pose const & pose, IndexList const & residues, ClosureProblemOP problem) {

	Real phi, psi;
	Ramachandran2B const & rama =
		ScoringManager::get_instance()->get_Ramachandran2B();

	BOOST_FOREACH ( Size residue, residues ) {

		// Currently we don't have data for both neighbors together.  Instead, for
		// each residue we consider either its left neighbor or its right neighbor
		// with equal probability.  The conditionals in this loop might be
		// problematic from a branch prediction perspective, because this is a
		// fairly inner loop, but this hasn't been actually benchmarked.

		if ( numeric::random::uniform() > 0.5 ) {
			rama.random_phipsi_from_rama_left(
				pose.aa(residue - 1), pose.aa(residue), phi, psi);
		} else {
			rama.random_phipsi_from_rama_right(
				pose.aa(residue), pose.aa(residue + 1), phi, psi);
		}

		if ( pose.residue( residue ).has_property( "D_AA" ) ) {
			phi *= -1.0;
			psi *= -1.0;
		}

		problem->perturb_phi(residue, phi, DEGREES);
		problem->perturb_psi(residue, psi, DEGREES);
	}
}

}
}
}
