// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/scoring/methods/InterchainEnvEnergy.cc
/// @brief  Statistically derived rotamer pair potentials
/// @details For docking (or between chains) only those residues at the interface
///      and between the two interfaces need to be evaluated
/// @author Monica Berrondo


// Unit headers
#include <protocols/scoring/methods/InterchainEnvEnergy.hh>
#include <protocols/scoring/methods/InterchainEnvEnergyCreator.hh>

// Package headers
#include <protocols/scoring/InterchainPotential.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnvPairPotential.hh>
#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/vector1.hh>

// Basic headers
#include <basic/Tracer.hh>

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace scoring {
namespace methods {

static THREAD_LOCAL basic::Tracer TR( "protocols.scoring.methods.InterchainEnvEnergy" );

/// @details This must return a fresh instance of the InterchainEnvEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
InterchainEnvEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return core::scoring::methods::EnergyMethodOP( new InterchainEnvEnergy );
}

core::scoring::ScoreTypes
InterchainEnvEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( interchain_env );
	sts.push_back( interchain_contact );
	return sts;
}


/// c-tor
InterchainEnvEnergy::InterchainEnvEnergy() :
	parent( core::scoring::methods::EnergyMethodCreatorOP( new InterchainEnvEnergyCreator ) ),
	interchain_potential_( * InterchainPotential::get_instance() ),
	env_potential_( core::scoring::ScoringManager::get_instance()->get_EnvPairPotential() )
{}


/// clone
core::scoring::methods::EnergyMethodOP
InterchainEnvEnergy::clone() const
{
	return core::scoring::methods::EnergyMethodOP( new InterchainEnvEnergy() );
}


void
InterchainEnvEnergy::setup_for_scoring( core::pose::Pose & pose, core::scoring::ScoreFunction const & ) const
{
	// compute interpolated number of neighbors at various distance cutoffs
	env_potential_.compute_centroid_environment( pose );
	interchain_potential_.compute_interface( pose );
	pose.update_residue_neighbors();
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
/// @brief calculate the env_score for residues at the interface
void
InterchainEnvEnergy::residue_energy(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose,
	core::scoring::EnergyMap & emap
) const
{
	using namespace core::scoring;
	core::Real env_score ( 0.0 );
	interchain_potential_.evaluate_env_score( pose, rsd, env_score );

	emap[ interchain_env ] += env_score;
}

// is there a better way to get the contact score??
void
InterchainEnvEnergy::finalize_total_energy(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & emap
) const
{
	using namespace core;
	using namespace core::scoring;
	Real contact_score ( 0.0 );
	//int interface_residues ( 0 );
	//interface_residues = potential_.interface_residues( pose );
	//TR.Debug << "residues at interface: " << interface_residues << std::endl;

	//contact_score = ( 20 - interface_residues ) * 0.5;
	//if ( interface_residues == 0 ) contact_score += 2.0;
	//if ( interface_residues == 1 ) contact_score += 1.0;
	//if ( interface_residues == 2 ) contact_score += 0.5;

	interchain_potential_.evaluate_contact_score( pose, contact_score );
	emap[ interchain_contact ] = contact_score;

	// sets calculated from the CenPairInfo to false
	interchain_potential_.finalize( pose );
}
core::Size
InterchainEnvEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
}
