// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/EnvPairEnergy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/util.hh>
#include <core/energy_methods/EnvEnergy.hh>
#include <core/energy_methods/EnvEnergyCreator.hh>

// Package headers
#include <core/scoring/EnvPairPotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/VariantType.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/EnergyMap.hh>


// Utility headers

// C++

namespace core {
namespace energy_methods {



/// @details This must return a fresh instance of the EnvEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
EnvEnergyCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const &
) const {
	return utility::pointer::make_shared< EnvEnergy >();
}

core::scoring::ScoreTypes
EnvEnergyCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( env );
	sts.push_back( cbeta );
	return sts;
}


/// c-tor
EnvEnergy::EnvEnergy() :
	parent( utility::pointer::make_shared< EnvEnergyCreator >() ),
	potential_( core::scoring::ScoringManager::get_instance()->get_EnvPairPotential() )
{}


/// clone
core::scoring::methods::EnergyMethodOP
EnvEnergy::clone() const
{
	return utility::pointer::make_shared< EnvEnergy >();
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
EnvEnergy::setup_for_scoring( pose::Pose & pose, core::scoring::ScoreFunction const & ) const
{
	// compute interpolated number of neighbors at various distance cutoffs
	pose.update_residue_neighbors();
	potential_.compute_centroid_environment( pose );
}


///////////////////////////////////////
//
// ENV SCORE AND CBETA SCORE
void
EnvEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	core::scoring::EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) return;
	if ( rsd.aa()==core::chemical::aa_unk ) return;

	Real env_score( 0.0 ), cb_score6( 0.0 ), cb_score12( 0.0 ), cb_score( 0.0 );

	potential_.evaluate_env_and_cbeta_scores( pose, rsd,
		env_score, cb_score6, cb_score12 );

	env_score *= 2.019;
	cb_score = 2.667 * ( cb_score6 + cb_score12 ) * 0.3;

	core::Real rsd_wt = core::scoring::methods::get_residue_weight_by_ss( pose.conformation().secstruct( rsd.seqpos() ) );

	emap[ core::scoring::env   ] += env_score * rsd_wt;
	emap[ core::scoring::cbeta ] += cb_score  * rsd_wt;
} // residue_energy

void
EnvEnergy::finalize_total_energy(
	pose::Pose & pose,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap &
) const
{
	potential_.finalize( pose );
}
core::Size
EnvEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
