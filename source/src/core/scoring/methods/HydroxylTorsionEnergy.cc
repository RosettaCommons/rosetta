// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/HydroxylTorsionEnergy.hh
/// @brief  Term for chi3 on tyrosine residues to prefer the hydrogen lie in the plane of the ring
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit headers
#include <core/scoring/methods/HydroxylTorsionEnergy.hh>
#include <core/scoring/methods/HydroxylTorsionEnergyCreator.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/constants.hh>
#include <numeric/deriv/dihedral_deriv.hh>
#include <core/scoring/DerivVectorPair.hh>

#include <basic/Tracer.hh>
#include <iostream>
#include <core/kinematics/Jump.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <basic/database/open.hh>

namespace core {
namespace scoring {
namespace methods {

static THREAD_LOCAL basic::Tracer TR( "core.scoring.HydroxylTorsionEnergy" );

/// @details This must return a fresh instance of the P_AA_pp_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
HydroxylTorsionEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new HydroxylTorsionEnergy );
}

ScoreTypes
HydroxylTorsionEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( hxl_tors );
	return sts;
}

/// ctor
HydroxylTorsionEnergy::HydroxylTorsionEnergy() :
	parent( methods::EnergyMethodCreatorOP( new HydroxylTorsionEnergyCreator ) )
{
	// hard-coded for now
	read_database( "scoring/score_functions/bondlength_bondangle/hydroxyl_torsion_preference.txt" );
}

/// clone
EnergyMethodOP
HydroxylTorsionEnergy::clone() const
{
	return EnergyMethodOP( new HydroxylTorsionEnergy );
}

// copied over from cart_bonded
std::string
HydroxylTorsionEnergy::get_restag( core::chemical::ResidueType const & restype ) const
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
HydroxylTorsionEnergy::read_database(
	std::string filename
)
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

void
HydroxylTorsionEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const
{
	std::string restag( get_restag( rsd.type() ) );

	tors_iterator it = torsion_params_.find( restag );
	if ( it == torsion_params_.end() ) return;

	Size const natm( rsd.natoms() );

	for ( it = torsion_params_.begin(); it != torsion_params_.end();
			it = torsion_params_.equal_range(it->first).second ) {

		if ( it->first != restag ) continue;

		if ( torsion_params_.count( it->first ) <= 1 ) continue;

		std::pair< tors_iterator, tors_iterator > range = torsion_params_.equal_range( restag );
		tors_iterator it2;

		for ( it2 = range.first; it2 != range.second; ++it2 ) {
			TorsionParams const &p = it2->second;

			debug_assert( (p.atm[1] <= natm) && (p.atm[2] <= natm) &&
				(p.atm[3] <= natm) && ( p.atm[4] <=natm ) );

			Real tors = numeric::dihedral_radians(
				rsd.xyz( p.atm[1] ), rsd.xyz( p.atm[2] ), rsd.xyz( p.atm[3] ), rsd.xyz( p.atm[4] ) );

			Real score = p.k * ( std::cos( p.n*tors - p.delta ) + 1 );
			emap[ hxl_tors ] += score;
			/*
			printf( "Atms/k/n/delta/tors/score: %2d %2d %2d %2d %5.2f %4.1f %6.4f %6.4f %8.3f\n",
			int(p.atm[1]), int(p.atm[2]), int(p.atm[3]), int(p.atm[4]),
			p.k, p.n, p.delta, tors, score
			);
			*/
		}
	}
}

// using atom derivs instead of dof derivates;
void
HydroxylTorsionEnergy::eval_residue_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & /*res_data_cache*/,
	pose::Pose const & /*pose*/,
	EnergyMap const & /*weights*/,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{
	std::string restag = get_restag( rsd.type() );

	Size const natm( rsd.natoms() );

	tors_iterator it = torsion_params_.find( restag );
	if ( it == torsion_params_.end() ) return;

	for ( it = torsion_params_.begin(); it != torsion_params_.end();
			it = torsion_params_.equal_range(it->first).second ) {

		if ( it->first != restag ) continue;
		if ( torsion_params_.count( it->first ) <= 1 ) continue;

		tors_iterator it2;
		std::pair< tors_iterator, tors_iterator > range = torsion_params_.equal_range( restag );

		for ( it2 = range.first; it2 != range.second; ++it2 ) {
			TorsionParams const &p = it2->second;

			debug_assert( (p.atm[1] <= natm) && (p.atm[2] <= natm) &&
				(p.atm[3] <= natm) && ( p.atm[4] <=natm ) );

			Real tors = numeric::dihedral_radians(
				rsd.xyz( p.atm[1] ), rsd.xyz( p.atm[2] ), rsd.xyz( p.atm[3] ), rsd.xyz( p.atm[4] ) );
			Real dE_dtors = -p.k * p.n *( std::sin( p.n*tors - p.delta ) );

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

	return;
}

void
HydroxylTorsionEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const {}

core::Size
HydroxylTorsionEnergy::version() const
{
	return 1; // Initial versioning
}

} // methods
} // scoring
} // core

