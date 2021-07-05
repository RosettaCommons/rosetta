// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/HydroxylTorsionEnergy.hh
/// @brief  Term for proton_chi on Ser/Thr/Tyr residues
/// @author Hahnbeom Park (hahnbeom@gmail.com)

// Unit headers
#include <core/scoring/HydroxylTorsionPotential.hh>
#include <core/energy_methods/HydroxylTorsionEnergy.hh>
#include <core/energy_methods/HydroxylTorsionEnergyCreator.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#include <basic/Tracer.hh>

namespace core {
namespace energy_methods {


static basic::Tracer TR( "core.energy_methods.HydroxylTorsionEnergy" );

/// @details This must return a fresh instance of the P_AA_pp_Energy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
HydroxylTorsionEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< HydroxylTorsionEnergy >();
}

core::scoring::ScoreTypes
HydroxylTorsionEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( hxl_tors );
	return sts;
}

/// ctor
HydroxylTorsionEnergy::HydroxylTorsionEnergy() :
	parent( utility::pointer::make_shared< HydroxylTorsionEnergyCreator >() ),
	potential_( core::scoring::ScoringManager::get_instance()->get_HydroxylTorsionPotential() )
{
	// hard-coded for now
}

/// clone
core::scoring::methods::EnergyMethodOP
HydroxylTorsionEnergy::clone() const
{
	return utility::pointer::make_shared< HydroxylTorsionEnergy >();
}

void
HydroxylTorsionEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	core::scoring::EnergyMap & emap
) const
{
	Real score = potential_.eval_residue_energy( rsd );
	emap[ core::scoring::hxl_tors ] += score;
}

// using atom derivs instead of dof derivates;
void
HydroxylTorsionEnergy::eval_residue_derivatives(
	conformation::Residue const & rsd,
	core::scoring::ResSingleMinimizationData const & /*res_data_cache*/,
	pose::Pose const & /*pose*/,
	core::scoring::EnergyMap const & /*weights*/,
	utility::vector1< core::scoring::DerivVectorPair > & atom_derivs
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

} // scoring
} // core

