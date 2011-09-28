// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RG_Energy_Fast.cc
/// @brief  Radius of gyration energy function definition.
/// @author James Thompson


// Unit headers
#include <core/scoring/methods/RG_Energy_Fast.hh>
#include <core/scoring/methods/RG_Energy_FastCreator.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

// AUTO-REMOVED #include <core/scoring/EnvPairPotential.hh>
//#include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>

// Project headers
#include <core/pose/Pose.hh>
// Auto-header: duplicate removed #include <core/conformation/Residue.hh>


// Utility headers
#include <basic/prof.hh>

//Auto Headers
#include <core/scoring/EnergyMap.hh>



// C++


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the RG_Energy_Fast class,
/// never an instance already in use
methods::EnergyMethodOP
RG_Energy_FastCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new RG_Energy_Fast;
}

ScoreTypes
RG_Energy_FastCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rg );
	return sts;
}


/// c-tor
RG_Energy_Fast::RG_Energy_Fast() :
	parent( new RG_Energy_FastCreator )
{}


/// clone
EnergyMethodOP
RG_Energy_Fast::clone() const
{
	return new RG_Energy_Fast;
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
RG_Energy_Fast::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	using namespace conformation;

	PROF_START( basic::RG );

	totals[ rg ] = calculate_rg_score( pose );

	PROF_STOP( basic::RG );
} // finalize_total_energy

core::Real
RG_Energy_Fast::calculate_rg_score( core::pose::Pose const & pose ) const
{

	Size const nres( pose.total_residue() );
	Size nres_counted=0;

	///////////////////////////////////////
	//
	// RG SCORE

	// calculate center of mass
	Vector center_of_mass( 0, 0, 0 );
	for ( Size i = 1; i <= nres; ++i ) {
		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( "REPLONLY" ) ) continue;
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;

		Vector const v( pose.residue(i).nbr_atom_xyz() );
		center_of_mass += v;
		nres_counted++;
	}
	center_of_mass /= nres_counted;

	// calculate RG based on distance from center of mass
	Real rg_score = 0;
	for ( Size i = 1; i <= nres; ++i ) {
		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( "REPLONLY" ) ) continue;
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;

		Vector const v( pose.residue(i).nbr_atom_xyz() );
		rg_score += v.distance_squared( center_of_mass );
	}

	// This definition of rg differs with the conventional definition which
	// divides by nres and not by nres-1.  For the sake of matching r++, it's
	// being left at nres-1 for now, but is a candidate for change in the near
	// future.
	rg_score /= (nres_counted - 1);

	return sqrt( rg_score );

}

core::Real
RG_Energy_Fast::calculate_rg_score(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & relevant_residues) const
{

	Size const nres( pose.total_residue() );
	Size nres_counted=0;

	///////////////////////////////////////
	//
	// RG SCORE

	// calculate center of mass
	Vector center_of_mass( 0, 0, 0 );
	for ( Size i = 1; i <= nres; ++i ) {
		if (!relevant_residues[i]) continue;

		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( "REPLONLY" ) ) continue;
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;

		Vector const v( pose.residue(i).nbr_atom_xyz() );
		center_of_mass += v;
		nres_counted++;
	}
	center_of_mass /= nres_counted;

	// calculate RG based on distance from center of mass
	Real rg_score = 0;
	for ( Size i = 1; i <= nres; ++i ) {
		if (!relevant_residues[i]) continue;

		// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
		if ( pose.residue(i).has_variant_type( "REPLONLY" ) ) continue;
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;

		Vector const v( pose.residue(i).nbr_atom_xyz() );
		rg_score += v.distance_squared( center_of_mass );
	}

	// This definition of rg differs with the conventional definition which
	// divides by nres and not by nres-1.  For the sake of matching r++, it's
	// being left at nres-1 for now, but is a candidate for change in the near
	// future.
	rg_score /= (nres_counted - 1);

	return sqrt( rg_score );

}

core::Size
RG_Energy_Fast::version() const
{
	return 1; // Initial versioning
}

} // methods
} // scoring
} // core
