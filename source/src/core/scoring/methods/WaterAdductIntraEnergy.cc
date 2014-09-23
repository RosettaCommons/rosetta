// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/WaterAdductIntraEnergy.hh
/// @brief  Energetic offset/cost for placing a water adduct on an amino or nucleic acid
/// @author Jim Havranek

// Unit headers
#include <core/scoring/methods/WaterAdductIntraEnergy.hh>
#include <core/scoring/methods/WaterAdductIntraEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyMap.hh>

// Project headers
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>



namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the WaterAdductIntraEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
WaterAdductIntraEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new WaterAdductIntraEnergy );
}

ScoreTypes
WaterAdductIntraEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( h2o_intra );
	return sts;
}


/// ctor
WaterAdductIntraEnergy::WaterAdductIntraEnergy() :
	parent( methods::EnergyMethodCreatorOP( new WaterAdductIntraEnergyCreator ) )
{}

/// clone
EnergyMethodOP
WaterAdductIntraEnergy::clone() const
{
	return EnergyMethodOP( new WaterAdductIntraEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

///
void
WaterAdductIntraEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const & ,
	EnergyMap & emap
) const
{
	// Sum over all waters
	for( int atm = 1, atme = rsd.natoms() ; atm <= atme ; ++atm ) {
		if( rsd.atom_type( atm ).is_h2o() ) emap[ h2o_intra ] += 1.0;
	}
}


///
Real
WaterAdductIntraEnergy::eval_dof_derivative(
	id::DOF_ID const &,// dof_id,
	id::TorsionID const &, //  tor_id
	pose::Pose const &, // pose
	ScoreFunction const &, //sfxn
	EnergyMap const & // weights
) const
{
	return 0.0;
}

/// @brief WaterAdductIntraEnergy is context independent; indicates that no
/// context graphs are required
void
WaterAdductIntraEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}
core::Size
WaterAdductIntraEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

