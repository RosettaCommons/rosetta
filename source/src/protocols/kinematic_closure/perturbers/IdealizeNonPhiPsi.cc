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
#include <protocols/kinematic_closure/internal.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/perturbers/IdealizeNonPhiPsi.hh>

// Core headers
#include <core/pose/Pose.hh>

// Utility headers
#include <numeric/conversions.hh>
#include <boost/foreach.hpp>

#define foreach BOOST_FOREACH

namespace protocols {
namespace kinematic_closure {
namespace perturbers {

using numeric::conversions::DEGREES;

void IdealizeNonPhiPsi::perturb_subset(
		Pose const &, IndexList const & residues, ClosureProblemOP problem) {

	foreach (Size residue, residues) {
		problem->perturb_omega(residue, IdealParameters::omega_dihedral, DEGREES);
		problem->perturb_c_n_ca(residue, IdealParameters::c_n_ca_angle, DEGREES);
		problem->perturb_n_ca_c(residue, IdealParameters::n_ca_c_angle, DEGREES);
		problem->perturb_ca_c_n(residue, IdealParameters::ca_c_n_angle, DEGREES);
		problem->perturb_c_n(residue, IdealParameters::c_n_length);
		problem->perturb_n_ca(residue, IdealParameters::n_ca_length);
		problem->perturb_ca_c(residue, IdealParameters::ca_c_length);
	}
}

void IdealizeNonPhiPsi::perturb_subset_with_balance(
		Pose const & pose, IndexList const & residues, ClosureProblemOP problem) {

	perturb_subset(pose, residues, problem);
}

}
}
}
