// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/FreeMoietyEnergy.cc
/// @brief  FreeMoiety energy method implementation
/// @author Rhiju Das (rhiju@stanford.edu)

// Unit headers
#include <core/scoring/methods/FreeMoietyEnergy.hh>
#include <core/scoring/methods/FreeMoietyEnergyCreator.hh>

// Package Headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>

// Project headers
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the FreeMoietyEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
FreeMoietyEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new FreeMoietyEnergy;
}

ScoreTypes
FreeMoietyEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( free_P );
	sts.push_back( free_2HOstar );
	return sts;
}



/// ctor
FreeMoietyEnergy::FreeMoietyEnergy() :
	parent( new FreeMoietyEnergyCreator )
{}

FreeMoietyEnergy::~FreeMoietyEnergy() {}

/// clone
core::scoring::methods::EnergyMethodOP
FreeMoietyEnergy::clone() const
{
	return new FreeMoietyEnergy;
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

/// @details Allocate the scratch space object on the stack to
/// alieviate thread-safety concerns.  Scratch does not use new.
void
FreeMoietyEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const
{

	if ( rsd.has_variant_type( "VIRTUAL_PHOSPHATE" ) )	      emap[ free_P ] += -1.0;
	if ( rsd.has_variant_type( "VIRTUAL_O2STAR_HYDROGEN" ) )	emap[ free_2HOstar ] += -1.0;

}


/// @brief FreeMoietyEnergy is context independent; indicates that no context graphs are required
void
FreeMoietyEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
) const
{}

core::Size
FreeMoietyEnergy::version() const
{
	return 1; // Initial versioning
}



} // methods
} // scoring
} // core

