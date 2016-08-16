// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/DirectReadoutEnergy.cc
/// @brief  Statistically derived DNA contact potential class implementation
/// @author Phil Bradley
/// @author Amy Ticoll

// C++ headers
#include <iostream>
#include <sstream>

// Unit headers
#include <core/scoring/methods/DirectReadoutEnergy.hh>
#include <core/scoring/methods/DirectReadoutEnergyCreator.hh>

// Package headers
#include <core/scoring/dna/DirectReadoutPotential.hh>
#include <core/scoring/ScoringManager.hh>
//#include <core/scoring/etable/count_pair/CountPair1BC4.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>


using namespace std;

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the DirectReadoutEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
DirectReadoutEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new DirectReadoutEnergy );
}

ScoreTypes
DirectReadoutEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( dna_dr );
	return sts;
}


/// @details  C-TOR
DirectReadoutEnergy::DirectReadoutEnergy() :
	parent( methods::EnergyMethodCreatorOP( new DirectReadoutEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_DirectReadoutPotential() )
{}


/// @details  Clone
EnergyMethodOP
DirectReadoutEnergy::clone() const
{
	return EnergyMethodOP( new DirectReadoutEnergy() );
}

/// @details  Totally inefficient implementation to avoid defining nbr-ness

void
DirectReadoutEnergy::finalize_total_energy(
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap
) const
{
	Size const nres( pose.total_residue() );

	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd1( pose.residue(i) );
		if ( !rsd1.is_protein() ) continue;

		for ( Size j=1; j<= nres; ++j ) {
			conformation::Residue const & rsd2( pose.residue(j) );

			if ( !rsd2.is_DNA() ) continue;

			my_residue_pair_energy( rsd1, rsd2, pose, scorefxn, emap );
		}
	}
}


/// @details  Unused rsd-pair implementation
void
DirectReadoutEnergy::my_residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const
{

	Real score(0.0);

	if ( rsd1.is_protein() && rsd2.is_DNA() ) {
		score = potential_.rsd_rsd_energy( rsd1, rsd2 );
	} else if ( rsd1.is_DNA() && rsd2.is_protein() ) {
		score = potential_.rsd_rsd_energy( rsd2, rsd1 );
	}

	emap[ dna_dr ] += score; // change

}
core::Size
DirectReadoutEnergy::version() const
{
	return 1; // Initial versioning
}


} // ns methods
} // ns scoring
} // ns core
