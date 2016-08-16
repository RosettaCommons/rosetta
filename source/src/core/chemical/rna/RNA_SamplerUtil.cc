// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/rna/RNA_SamplerUtil.cc
/// @author Rhiju Das

// Unit headers
#include <core/chemical/rna/RNA_SamplerUtil.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/types.hh>

// Package headers
#include <utility/vector1.hh>

// Utility headers
#include <numeric/angle.functions.hh>

// C++

namespace core {
namespace chemical {
namespace rna {

/////////////////////////////////////////////////////////////////////////////
// used in screening and packing
void add_values_from_center(
	utility::vector1<core::Real> & torsions,
	Real const center,
	Real const max_range,
	Real const bin_size
) {
	Real const epsilon = 0.0001; //A small number to avoid numerical error for 'Real'.
	for ( Real offset = 0; offset <= max_range + epsilon; offset += bin_size ) {
		torsions.push_back( center + offset );
		if ( offset != 0 ) torsions.push_back( center - offset );
	}
	std::sort( torsions.begin(), torsions.end() );
}

/////////////////////////////////////////////////////////////////////////////
// used in screening & packing
utility::vector1< Real >
get_full_torsions( Real const bin_size /* = 20.0 */ ){
	utility::vector1< Real > full_torsions;
	add_values_from_center( full_torsions, 0, 180, bin_size );
	//Avoid sampling both -180 and 180 deg
	Real const epsil =  0.001; //Arbitary small number
	if ( 360 - ( full_torsions.back() - full_torsions.front() ) < epsil ) full_torsions.pop_back();
	// TR << "FULL TORSIONS: "  << full_torsions << std::endl;
	return full_torsions;
}

/////////////////////////////////////////////////////////////////////////////
utility::vector1< Real >
get_epsilon_torsions( Real const delta,
	bool const extra_epsilon,
	Real const bin_size /* = 20.0*/ ){
	static chemical::rna::RNA_FittedTorsionInfo const torsion_info;
	return get_epsilon_torsions( numeric::principal_angle_degrees( delta ) <= torsion_info.delta_cutoff(),
		extra_epsilon, bin_size );
}

/////////////////////////////////////////////////////////////////////////////
utility::vector1< Real >
get_epsilon_torsions( bool const north_pucker,
	bool const extra_epsilon,
	Real const /* bin_size=20.0 */ )
{
	using namespace chemical::rna;
	static chemical::rna::RNA_FittedTorsionInfo const torsion_info;

	/////Epsilon rotamers/////
	//default: center +- 20 deg
	//extra_epsilon: center +- 60 deg
	utility::vector1< Real > epsilon_torsions;
	Real center = ( north_pucker ) ? torsion_info.epsilon_north() : torsion_info.epsilon_south();
	Real max_range = 20.0;
	if ( extra_epsilon ) {
		max_range = 60.0;
		//Choice made by Parin, to cover the uneven tails of
		//epsilons distributions by extra_epsilon mode.
		if ( north_pucker ) {
			center -= 20.0;
		} else {
			center += 20.0;
		}
	}
	Real const bin_size_( 10.0 );
	add_values_from_center( epsilon_torsions, center, max_range, bin_size_ );
	// TR << "EPSILON CENTER: "  << center << "   based on DELTA: " << delta << std::endl;
	// TR << "EPSILON TORSIONS: "  << epsilon_torsions << std::endl;
	return epsilon_torsions;
}

} //ns rna
} //ns chemical
} //ns core
