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
#include <protocols/kinematic_closure/perturbers/OmegaPerturber.hh>

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

void OmegaPerturber::perturb_subset(
	Pose const & pose, IndexList const & residues, ClosureProblemOP problem) {

	using core::chemical::aa_pro;
	using numeric::random::gaussian;
	using numeric::random::uniform;
	using numeric::conversions::DEGREES;

	BOOST_FOREACH ( Size residue, residues ) {
		// Omega distribution mean and stddev from Berkholz et al., PNAS 2012.
		Real trans_omega = 179.1 + 6.3 * gaussian();

		// There's very little variation, currently not captured here at all.
		Real cis_omega = 0;

		// Pick which omega to apply.  For pre-prolines, the cis torsion is picked
		// 0.1% of the time.  For everything else , the trans is always picked.
		Real cis_prob = (pose.aa(residue + 1) == aa_pro) ? 0.001 : 0;
		Real omega = (cis_prob > uniform()) ? cis_omega : trans_omega;

		// Apply the omega to the closure problem.
		problem->perturb_omega(residue, omega, DEGREES);
	}
}

}
}
}
