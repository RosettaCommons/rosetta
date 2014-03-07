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
#include <protocols/kinematic_closure/perturbers/RamaPerturber.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/Ramachandran.hh>
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

using core::scoring::Ramachandran;
using core::scoring::ScoringManager;
using numeric::conversions::DEGREES;

void RamaPerturber::perturb_subset(
		Pose const & pose, IndexList const & residues, ClosureProblemOP problem) {

	Real phi, psi;
	Ramachandran const & rama =
		ScoringManager::get_instance()->get_Ramachandran();

	BOOST_FOREACH(Size residue, residues) {
		rama.random_phipsi_from_rama(pose.aa(residue), phi, psi);
		problem->perturb_phi(residue, phi, DEGREES);
		problem->perturb_psi(residue, psi, DEGREES);
	}
}

/// @details The balanced version of this algorithm picks phi/psi pairs 
/// uniformly from the allowed regions of rama space.  This is more efficient 
/// than sampling from a completely uniform distribution, and it still obeys 
/// detailed balance.

void RamaPerturber::perturb_subset_with_balance(
		Pose const & pose, IndexList const & residues, ClosureProblemOP problem) {

	Real phi, psi;
	Ramachandran const & rama =
		ScoringManager::get_instance()->get_Ramachandran();

	BOOST_FOREACH(Size residue, residues) {
		rama.uniform_phipsi_from_allowed_rama(pose.aa(residue), phi, psi);
		problem->perturb_phi(residue, phi, DEGREES);
		problem->perturb_psi(residue, psi, DEGREES);
	}
}

}
}
}
