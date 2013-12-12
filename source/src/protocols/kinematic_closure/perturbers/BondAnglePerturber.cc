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
#include <protocols/kinematic_closure/perturbers/BondAnglePerturber.hh>

// Core headers
#include <core/pose/Pose.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/kinematic_closure/vector.hh>
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

namespace protocols {
namespace kinematic_closure {
namespace perturbers {

void BondAnglePerturber::perturb_subset(
		Pose const & pose, IndexList const & residues, ClosureProblemOP problem) {

	using numeric::random::gaussian;
	using numeric::conversions::DEGREES;

	// Sorry for using magic numbers.  Both values were taken from a gaussian fit 
	// of the bond angle distribution observed in the Top8000 database.

	foreach(Size residue, residues) {
		Real angle = 111.24096 + 2.27632 * gaussian();
		problem->perturb_n_ca_c(residue, angle, DEGREES);
	}
}

}
}
}
