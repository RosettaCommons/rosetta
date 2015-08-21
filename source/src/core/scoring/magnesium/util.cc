// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/magnesium/util.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/magnesium/util.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "core.scoring.magnesium.util" );

using namespace core::chemical::rna;

namespace core {
namespace scoring {
namespace magnesium {

////////////////////////////////////////////////////////////////////////////////////////////////////////
// will deprecate this one soon...
Real
get_cos_theta( core::conformation::Residue const & rsd1,
	Size const i,  Vector const & j_xyz,
	Size const i_base /* = 0 */ )
{
	Size i_base_local( i_base );
	Vector base_xyz( 0.0 );
	return get_cos_theta( rsd1, i, j_xyz, i_base_local, base_xyz );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
get_cos_theta( core::conformation::Residue const & rsd1,
	Size const i,  Vector const & j_xyz,
	Size & i_base,
	Vector & xyz_base )
{

	static hbonds::HBondOptions const hbond_options;

	Vector const & i_xyz = rsd1.xyz( i );

	Vector dummy;
	if ( i_base > 0 ) {
		xyz_base = rsd1.xyz( i_base );
	} else {
		i_base = rsd1.atom_base( i );
		if ( rsd1.heavyatom_is_an_acceptor( i ) ) {
			chemical::Hybridization acc_hybrid( rsd1.atom_type( i ).hybridization() );
			make_hbBasetoAcc_unitvector(
				hbond_options,
				acc_hybrid,
				rsd1.atom( i ).xyz(),
				rsd1.xyz( i_base ),
				rsd1.xyz( rsd1.abase2( i ) ),
				xyz_base, dummy );
		} else if ( rsd1.atom_type( i ).name() == "Hpol" ) {
			xyz_base = rsd1.xyz( i_base );
		} else {
			return -999; // bogus value
		}
	}

	Vector const a = j_xyz - i_xyz;
	Vector const b = xyz_base - i_xyz;

	Real cos_theta = dot( a, b )/( a.length() * b.length() );

	//  if ( a.length() < 3.0 ) std::cout << rsd1.name1() << " " << rsd1.seqpos() << " " << rsd1.atom_name(i) << " " << cos_theta << std::endl;

	return cos_theta;
}



/////////////////////////////////
Real
get_gaussian_potential_score(
	GaussianParameter const & mg_potential_gaussian_parameter,
	Vector const & pos1,
	Vector const & pos2 )
{ // later expand to do derivative calculation

	Distance const d = ( pos1 - pos2 ).length();

	return get_gaussian_score( mg_potential_gaussian_parameter, d );
}

///////////////////////////////////////////////////
Real
get_gaussian_score(
	GaussianParameter const & mg_potential_gaussian_parameter,
	Real const & d )
{ // later expand to do derivative calculation

	Real const a     = mg_potential_gaussian_parameter.amplitude;
	Real const d0    = mg_potential_gaussian_parameter.center;
	Real const sigma = mg_potential_gaussian_parameter.width;

	Real const score = a * exp( -0.5 * std::pow( ( d - d0 )/sigma, 2 ) );
	return score;
}


///////////////////////////////////////////////////
Real
get_gaussian_deriv(
	GaussianParameter const & mg_potential_gaussian_parameter,
	Real const & d )
{

	Real const a     = mg_potential_gaussian_parameter.amplitude;
	Real const d0    = mg_potential_gaussian_parameter.center;
	Real const sigma = mg_potential_gaussian_parameter.width;

	Real const deriv = a * ( d - d0 )/( sigma * sigma ) * exp( -0.5 * std::pow( ( d - d0 )/sigma, 2 ) );
	return deriv;
}


///////////////////////////////////////////////////
void
get_closest_orbital_axis( core::conformation::Residue const & mg_rsd, Vector const & j_xyz,
	Size & which_orbital_atom, Real & cos_angle ) {
	runtime_assert( mg_rsd.natoms() == 7 ); // sanity check that this is MG, V1, V2, ... V6
	Vector const & mg_xyz = mg_rsd.xyz( 1 );
	Vector lig = ( j_xyz - mg_xyz ).normalized();
	Real best_cos_angle( -1.0 );
	for ( Size n = 1; n <= 6; n++ ) {
		Vector v = ( mg_rsd.xyz( n+1 ) - mg_xyz ).normalized();
		Real const cos_angle = dot( v, lig ); // -1.0 (worst) to 1.0 (best)
		if ( cos_angle >= best_cos_angle ) {
			best_cos_angle = cos_angle;
			which_orbital_atom = n;
		}
	}
	cos_angle = best_cos_angle;
}

///////////////////////////////////////////////////
Real
get_cos_angle_to_closest_orbital_axis( core::conformation::Residue const & mg_rsd, Vector const & j_xyz ) {
	Size which_orbital_atom;
	Real cos_angle;
	get_closest_orbital_axis( mg_rsd, j_xyz, which_orbital_atom, cos_angle );
	return cos_angle;
}

///////////////////////////////////////////////////
Size
get_closest_orbital_axis( core::conformation::Residue const & mg_rsd, Vector const & j_xyz ) {
	Size which_orbital_atom( 0 );
	Real cos_angle( 0.0 );
	get_closest_orbital_axis( mg_rsd, j_xyz, which_orbital_atom, cos_angle );
	runtime_assert( which_orbital_atom >= 1 && which_orbital_atom <= 6 );
	return which_orbital_atom;
}



} //magnesium
} //scoring
} //core
