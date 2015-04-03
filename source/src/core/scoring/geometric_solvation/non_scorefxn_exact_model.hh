// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/geometric_solvation/exact_model.hh
/// @brief  Solvation model based on penalizing potential for Hbonding to solvent
/// @author John Karanicolas


#ifndef INCLUDED_core_scoring_geometric_solvation_non_scorefxn_exact_model_hh
#define INCLUDED_core_scoring_geometric_solvation_non_scorefxn_exact_model_hh

#include <core/types.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/geometric_solvation/ExactOccludedHbondSolEnergy.hh>

// Utility headers

// C++ headers
#include <map>

namespace core {
namespace scoring {
namespace geometric_solvation {

void add_to_individual_sol_energies( pose::Pose & input_pose,
	core::Size const polar_resnum,
	core::Size const polar_atomno,
	core::scoring::etable::EtableOP etable_ptr,
	GridInfo const & grid_info,
	core::Real const & grid_constant,
	std::vector < std::vector < std::vector <core::Real> > > const & water_weights,
	std::vector < std::vector < std::vector <bool> > > & occluded_sites,
	bool const hydrogens_can_occlude,
	bool const pairwise_additive,
	bool const pairwise_additive_output,
	utility::vector1 <core::Real> & residue_energies
);


core::Real compute_exact_geosol( pose::Pose & input_pose,
	bool const hydrogens_can_occlude,
	bool const pairwise_additive,
	bool const pairwise_additive_output,
	utility::vector1<core::Real> & residue_energies
);


} // geometric_solvation
} // scoring
} // core

#endif

