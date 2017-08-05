// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/HydroxylTorsionPotential.cc
/// @brief  Potential for core/scoring/methods/HydroxyTorsionEnergy
/// @author Hahnbeom Park (hahnbeom@gmail.com)

#include <core/scoring/HydroxylTorsionPotential.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/constants.hh>
#include <numeric/deriv/dihedral_deriv.hh>
#include <core/scoring/DerivVectorPair.hh>

#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <basic/database/open.hh>
#include <iostream>

namespace core {
namespace scoring {

static THREAD_LOCAL basic::Tracer TR("core.scoring.HydroxylTorsionPotential");

/// ctor
HydroxylTorsionPotential::HydroxylTorsionPotential()
{
	// hard-code for now....
	read_database( "scoring/score_functions/bondlength_bondangle/hydroxyl_torsion_preference.txt" );
}

// copied over from cart_bonded
std::string
HydroxylTorsionPotential::get_restag( core::chemical::ResidueType const & restype ) const
{
	using namespace core::chemical;

	if ( core::chemical::is_canonical_D_aa(restype.aa()) ) {
		std::string rsdname = restype.name();
		if ( rsdname.substr(0, rsdname.find(chemical::PATCH_LINKER)) == "DHIS_D" ) rsdname="HIS_D"; //If this is a DHIS_D, return HIS_D.
		else rsdname=core::chemical::name_from_aa( core::chemical::get_L_equivalent( restype.aa() ) ); //Otherwise, for D-amino acids, return the L-equivalent.
		return rsdname;
	} else if ( !restype.is_protein() && !restype.is_NA() ) {
		return restype.name3();
	} else {
		std::string rsdname = restype.name();
		rsdname = rsdname.substr( 0, rsdname.find(chemical::PATCH_LINKER) );
		return rsdname;
	}
}

void
HydroxylTorsionPotential::read_database( std::string filename )
{
	utility::io::izstream instream;
	basic::database::open( instream, filename );

	std::string fileline;
	std::string restag;
	std::string atm1, atm2, atm3, atm4;
	Real k, n, delta;

	chemical::ResidueTypeSetCOP rsd_set
		= chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	TR.Debug << "Reading database from " << filename << std::endl;

	while ( instream ) {
		getline(instream, fileline);
		std::istringstream linestream(fileline);

		linestream >> restag;
		if ( !rsd_set->has_name3( restag ) ) continue;

		linestream >> atm1 >> atm2 >> atm3 >> atm4
			>> k >> n >> delta;

		/*
		boost::unordered_multimap< std::string, TorsionParams >::const_iterator it
		= torsion_params_.find( restag );
		if( it == torsion_params_.end() ){
		torsion_params_[res3] = utility::vector1< TorsionParams >();
		torsion_params_.insert( std::make_pair( res3, utility::vector1< TorsionParams >() ) );
		TR.Debug << "add " << restag << std::endl;
		}
		*/

		chemical::ResidueType const &rsd = rsd_set->name_map( restag );

		TorsionParams p; p.atm.resize(4);
		p.atm[1] = rsd.atom_index( atm1 ); p.atm[2] = rsd.atom_index( atm2 );
		p.atm[3] = rsd.atom_index( atm3 ); p.atm[4] = rsd.atom_index( atm4 );
		p.k = k; p.n = n; p.delta = delta;
		torsion_params_.insert( std::make_pair( restag, p ) );
		TR.Debug << "added " << restag << ", total " << torsion_params_.size() << std::endl;
	}
}

Real
HydroxylTorsionPotential::eval_residue_energy(
	conformation::Residue const & rsd
) const
{
	Real score( 0.0 );

	bool const is_d( rsd.type().is_d_aa() );

	std::string restag = get_restag( rsd.type() );

	tors_iterator it = torsion_params_.find( restag );
	if ( it == torsion_params_.end() ) return 0.0;

	Size const natm( rsd.natoms() );

	for ( it = torsion_params_.begin(); it != torsion_params_.end();
			it = torsion_params_.equal_range(it->first).second ) {

		if ( it->first != restag ) continue;

		if ( torsion_params_.count( it->first ) == 0 ) continue;

		std::pair< tors_iterator, tors_iterator > range = torsion_params_.equal_range( restag );
		tors_iterator it2;

		for ( it2 = range.first; it2 != range.second; ++it2 ) {
			TorsionParams const &p = it2->second;

			debug_assert( (p.atm[1] <= natm) && (p.atm[2] <= natm) &&
				(p.atm[3] <= natm) && ( p.atm[4] <=natm ) );

			Real tors = numeric::dihedral_radians(
				rsd.xyz( p.atm[1] ), rsd.xyz( p.atm[2] ), rsd.xyz( p.atm[3] ), rsd.xyz( p.atm[4] ) );

			if ( is_d ) tors *= -1.0;

			score += p.k * ( std::cos( p.n*tors - p.delta ) + 1 );
			/*
			printf( "Atms/k/n/delta/tors/score: %2d %2d %2d %2d %5.2f %4.1f %6.4f %6.4f %8.3f\n",
			int(p.atm[1]), int(p.atm[2]), int(p.atm[3]), int(p.atm[4]),
			p.k, p.n, p.delta, tors, score
			);
			*/
		}
	}

	return score;
}

void
HydroxylTorsionPotential::eval_residue_derivative(
	conformation::Residue const & rsd,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{
	std::string restag = get_restag( rsd.type() );

	bool const is_d( rsd.type().is_d_aa() );

	Size const natm( rsd.natoms() );

	tors_iterator it = torsion_params_.find( restag );
	if ( it == torsion_params_.end() ) return;

	for ( it = torsion_params_.begin(); it != torsion_params_.end();
			it = torsion_params_.equal_range(it->first).second ) {

		if ( it->first != restag ) continue;
		if ( torsion_params_.count( it->first ) == 0 ) continue;

		tors_iterator it2;
		std::pair< tors_iterator, tors_iterator > range = torsion_params_.equal_range( restag );

		for ( it2 = range.first; it2 != range.second; ++it2 ) {
			TorsionParams const &p = it2->second;

			debug_assert( (p.atm[1] <= natm) && (p.atm[2] <= natm) &&
				(p.atm[3] <= natm) && ( p.atm[4] <=natm ) );

			Real tors = numeric::dihedral_radians(
				rsd.xyz( p.atm[1] ), rsd.xyz( p.atm[2] ), rsd.xyz( p.atm[3] ), rsd.xyz( p.atm[4] ) );
			if ( is_d ) tors *= -1.0;
			Real dE_dtors = -p.k * p.n *( std::sin( p.n*tors - p.delta ) );
			if ( is_d ) dE_dtors *= -1.0;

			Vector const &xyz1( rsd.xyz( p.atm[1] ) );
			Vector const &xyz2( rsd.xyz( p.atm[2] ) );
			Vector const &xyz3( rsd.xyz( p.atm[3] ) );
			Vector const &xyz4( rsd.xyz( p.atm[4] ) );

			Vector f1( 0.0 ), f2( 0.0 );
			numeric::deriv::dihedral_p1_cosine_deriv( xyz1, xyz2, xyz3, xyz4, tors, f1, f2 );
			atom_derivs[ p.atm[1] ].f1() += dE_dtors * f1;
			atom_derivs[ p.atm[1] ].f2() += dE_dtors * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv( xyz1, xyz2, xyz3, xyz4, tors, f1, f2 );
			atom_derivs[ p.atm[2] ].f1() += dE_dtors * f1;
			atom_derivs[ p.atm[2] ].f2() += dE_dtors * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p2_cosine_deriv( xyz4, xyz3, xyz2, xyz1, tors, f1, f2 );
			atom_derivs[ p.atm[3] ].f1() += dE_dtors * f1;
			atom_derivs[ p.atm[3] ].f2() += dE_dtors * f2;

			f1 = f2 = Vector(0.0);
			numeric::deriv::dihedral_p1_cosine_deriv( xyz4, xyz3, xyz2, xyz1, tors, f1, f2 );
			atom_derivs[ p.atm[4] ].f1() += dE_dtors * f1;
			atom_derivs[ p.atm[4] ].f2() += dE_dtors * f2;
		}
	}
}

} // namespace scoring
} // namespace core
