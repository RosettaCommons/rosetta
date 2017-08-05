// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/WaterSpecificEnergy.hh
/// @brief  Energetic considerations for explicit water molecules,
//
/// @author Joaquin Ambia, Jason K. Lai

// Unit headers
#include <core/scoring/methods/WaterSpecificEnergy.hh>
#include <core/scoring/methods/WaterSpecificEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.hh>

// Options header
#include <basic/options/option.hh>
#include <basic/options/keys/hydrate.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the WaterSpecificEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
WaterSpecificEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new WaterSpecificEnergy );
}

ScoreTypes
WaterSpecificEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( wat_desolv );
	return sts;
}


/// ctor
WaterSpecificEnergy::WaterSpecificEnergy() :
	parent( methods::EnergyMethodCreatorOP( new WaterSpecificEnergyCreator ) )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::hydrate;
	using namespace basic::options::OptionKeys::score;

	// hydrate/SPaDES protocol
	wat_desolv_ = option[ water_desolvation ]();
	if ( !option[ water_desolvation ].user() ) {
		if ( option[ water_hybrid_sf ]() ) {
			wat_desolv_ = 4.8;
		}
		if ( option[ short_residence_time_mode ] ) {
			wat_desolv_ = wat_desolv_ / 2.0;
		}
	}
}

/// clone
EnergyMethodOP
WaterSpecificEnergy::clone() const {
	return EnergyMethodOP( new WaterSpecificEnergy( *this ) );
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

///
void
WaterSpecificEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const
{
	// These terms are only meant to be used with water
	if ( rsd.name() != "TP3" ) return;

	// The term wat_desolv is the desolvation energy of a water, it is 0 if the molecules is in bulk water (away rotamers)
	// and it has a constant value if the molecule is near the protein
	// This might need to be revisited, for now we find away rotamers simply by their position,
	// We might actually want to check for neighbors
	if ( rsd.name() != "TP3" || rsd.xyz(1).x() > pose.residue(1).xyz(1).x() + 10000 ) emap[ wat_desolv ] = 0;
	else  emap[ wat_desolv ] = wat_desolv_;

}


///
Real
WaterSpecificEnergy::eval_dof_derivative(
	id::DOF_ID const &,// dof_id,
	id::TorsionID const &, //  tor_id
	pose::Pose const &, // pose
	ScoreFunction const &, //sfxn
	EnergyMap const & // weights
) const
{
	return 0.0;
}

/// @brief WaterSpecificEnergy is context independent; indicates that no
/// context graphs are required
void
WaterSpecificEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}
core::Size
WaterSpecificEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

