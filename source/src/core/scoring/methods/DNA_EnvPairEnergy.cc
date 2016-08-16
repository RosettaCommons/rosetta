// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/scoring/methods/DNA_EnvPairEnergy.cc
/// @brief  dna scoring
/// @author Phil Bradley
/// @author Amy Ticoll

// C++ headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>

// Unit headers
#include <core/scoring/methods/DNA_EnvPairEnergy.hh>
#include <core/scoring/methods/DNA_EnvPairEnergyCreator.hh>

// Package headers
#include <core/scoring/dna/DirectReadoutPotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <ObjexxFCL/format.hh>
#include <core/conformation/Residue.hh>


//using namespace std;

namespace core {
namespace scoring {
namespace methods {
methods::EnergyMethodOP
DNA_EnvPairEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & //options
) const {
	return EnergyMethodOP( new DNA_EnvPairEnergy() );
}

ScoreTypes
DNA_EnvPairEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( dna_env );
	sts.push_back( dna_pair );
	return sts;
}


/// VERSION 1
Size
DNA_EnvPairEnergy::version() const { return 1; }


/// @details  C-TOR
DNA_EnvPairEnergy::DNA_EnvPairEnergy(): // constructor
	parent( EnergyMethodCreatorOP( new DNA_EnvPairEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_DNA_EnvPairPotential() )
{
}


/// @details  Clone
EnergyMethodOP
DNA_EnvPairEnergy::clone() const
{
	return EnergyMethodOP( new DNA_EnvPairEnergy() );
}

/// @details  Totally inefficient implementation to avoid defining nbr-ness

void
DNA_EnvPairEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &, //scorefxn,
	EnergyMap & emap
) const
{
	//PROF_START( basic::DNA_ENV_PAIR_ENERGY );

	Size const nres( pose.total_residue() );
	Real const nbr_dis2_threshold( potential_.nbr_dis2_threshold() );

	/// get "centroid" positions
	utility::vector1< Vector > centroid_xyz( nres );
	for ( Size i=1; i<= nres; ++i ) centroid_xyz[ i ] = potential_.centroid_xyz( pose.residue(i) );

	///
	Real env_score( 0.0 ), pair_score( 0.0 );

	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd1( pose.residue(i) );
		chemical::AA const & aa( rsd1.aa() );

		if ( !rsd1.is_protein() ) continue;

		Size nbr_count(0);
		Vector const & xyz1( centroid_xyz[ i ] );

		for ( Size j=1; j<= nres; ++j ) {
			conformation::Residue const & rsd2( pose.residue(j) );

			if ( !rsd2.is_DNA() ) continue;

			Real const dis2( xyz1.distance_squared( centroid_xyz[ j ] ) );

			if ( dis2 < nbr_dis2_threshold ) {
				++nbr_count;
				pair_score += potential_.residue_pair_score( rsd2.aa(), aa, dis2 );
			}
		}
		env_score += potential_.residue_env_score( aa, nbr_count );
	}

	emap[ dna_env  ] = env_score;
	emap[ dna_pair ] = pair_score;

	//PROF_STOP( basic::DNA_ENV_PAIR_ENERGY );
}


} // ns methods
} // ns scoring
} // ns core
