// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/EnvPairEnergy.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/MembraneEnvEnergy.hh>
#include <core/scoring/methods/MembraneEnvEnergyCreator.hh>
// Package headers
#include <core/scoring/MembranePotential.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>


// Utility headers

// C++

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the MembraneEnvEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
MembraneEnvEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MembraneEnvEnergy );
}

ScoreTypes
MembraneEnvEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( Menv );
	return sts;
}


/// c-tor
MembraneEnvEnergy::MembraneEnvEnergy() :
	parent( methods::EnergyMethodCreatorOP( new MembraneEnvEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_MembranePotential() )
{}


/// clone
EnergyMethodOP
MembraneEnvEnergy::clone() const
{
	return EnergyMethodOP( new MembraneEnvEnergy );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////


void
MembraneEnvEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	// compute interpolated number of neighbors at various distance cutoffs
	pose.update_residue_neighbors();
	potential_.compute_centroid_environment( pose );
	potential_.compute_membrane_embedding( pose );

}

void
MembraneEnvEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sf) const
{
	setup_for_scoring(pose,sf);
}


///////////////////////////////////////
//
// ENV SCORE
void
MembraneEnvEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	Real env_score( 0.0 );
	// std::cout << rsd.seqpos() << "\n "; << depth << "\n";
	//std::cout << pos << "\n";
	if ( rsd.aa() == core::chemical::aa_vrt ) return;
	if ( MembraneTopology_from_pose( pose ).allow_scoring(rsd.seqpos()) ) {
		potential_.evaluate_env( pose, rsd, env_score);
	}

	emap[ Menv ] += env_score;
} // residue_energy

void
MembraneEnvEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap &
) const
{
	potential_.finalize( pose );
}
MembraneTopology const &
MembraneEnvEnergy::MembraneTopology_from_pose( pose::Pose const & pose ) const
{
	return *( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ));
}
core::Size
MembraneEnvEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
}
