// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/HydroxylTorsionEnergy.hh
/// @brief  Term for proton_chi on Ser/Thr/Tyr residues
/// @author Hahnbeom Park (hahnbeom@gmail.com)

// Unit headers
#include <core/scoring/HydroxylTorsionPotential.hh>
#include <core/scoring/methods/HydroxylTorsionEnergy.hh>
#include <core/scoring/methods/HydroxylTorsionEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <basic/Tracer.hh>

namespace core {
namespace scoring {
namespace methods {

static basic::Tracer TR( "core.scoring.HydroxylTorsionEnergy" );

/// @details This must return a fresh instance of the P_AA_pp_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
HydroxylTorsionEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new HydroxylTorsionEnergy );
}

ScoreTypes
HydroxylTorsionEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( hxl_tors );
	return sts;
}

/// ctor
HydroxylTorsionEnergy::HydroxylTorsionEnergy() :
	parent( methods::EnergyMethodCreatorOP( new HydroxylTorsionEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_HydroxylTorsionPotential() )
{
	// hard-coded for now
}

/// clone
EnergyMethodOP
HydroxylTorsionEnergy::clone() const
{
	return EnergyMethodOP( new HydroxylTorsionEnergy );
}

void
HydroxylTorsionEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const
{
	Real score = potential_.eval_residue_energy( rsd );
	emap[ hxl_tors ] += score;
}

// using atom derivs instead of dof derivates;
void
HydroxylTorsionEnergy::eval_residue_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & /*res_data_cache*/,
	pose::Pose const & /*pose*/,
	EnergyMap const & /*weights*/,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{
	potential_.eval_residue_derivative( rsd, atom_derivs );
	return;
}

void
HydroxylTorsionEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const {}

core::Size
HydroxylTorsionEnergy::version() const
{
	return 1; // Initial versioning
}

} // methods
} // scoring
} // core

