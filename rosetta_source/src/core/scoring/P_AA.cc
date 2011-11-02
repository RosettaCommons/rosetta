// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/P_AA.cc
/// @brief  Amino acid probability arrays and functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Andrew Leaver-Fay -- porting Stuarts code

// Unit headers
#include <core/scoring/P_AA.hh>

// Project headers
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>

// Numeric headers
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/interpolation/periodic_range/full/interpolation.hh>
#include <numeric/interpolation/periodic_range/half/interpolation.hh>

// ObjexxFCL headers
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

// Utility headers
#include <utility/io/izstream.hh>

// C++ headers
#include <cassert>

#include <utility/vector1.hh>



namespace core {
namespace scoring {


/// @brief Amino acid probability array: P(aa)
//Probability_AA P_AA;

/// @brief Amino acid conditional probability wrt number of neighbors array: P(aa|neighbors)
//Probability_AA_n P_AA_n;

/// @brief Amino acid conditional probability wrt (phi,psi) array: P(aa|phi,psi)
//Probability_AA_pp P_AA_pp;


/// @brief ctor -- Initialize the amino acid probability data structures
P_AA::P_AA()
{
	read_P_AA();
	read_P_AA_n();
	read_P_AA_pp();
}

P_AA::~P_AA() {}


/// @brief Read the amino acid probability file into P_AA
///
/// @note  Only the keys present in the file are given entries
void
P_AA::read_P_AA()
{
	using namespace core::chemical;

	// Read the probability file and load the array
	std::string id;
	Probability probability, probability_sum( 0.0 );
	utility::io::izstream stream;
	basic::database::open( stream, "P_AA" );

	P_AA_.resize( num_canonical_aas );

	while ( stream ) {
		using namespace ObjexxFCL::fmt;
		stream >> bite( 3, id ) >> skip( 1 ) >> bite( 9, probability ) >> skip;
		if ( stream ) {
			assert( ( probability >= Probability( 0.0 ) ) && ( probability <= Probability( 1.0 ) ) );
			AA aa = aa_from_name( id );
			assert( ( aa >= 1 ) && ( aa <= num_canonical_aas ) );
			probability_sum += probability;
			P_AA_[ aa ] = probability;
		} //! ADD INPUT ERROR HANDLING
	}
	stream.close();

	// Check probabilities sum to ~ 1
	assert( numeric::eq_tol( probability_sum, Probability( 1.0 ), Probability( .0001 ), Probability( .0001 ) ) );
}


/// @brief Read the amino acid conditional probability wrt (neighbors) file into P_AA_n
///
/// @note  Only the keys present in the file are given entries
/// @note  The file entries can be in any order
void
P_AA::read_P_AA_n()
{
	using namespace core::chemical;

	// Read the probability file and load the array
	std::string id;
	int n; // Number of neighbors
	Probability probability;
	utility::io::izstream stream;
	basic::database::open( stream, "P_AA_n" );

	P_AA_n_.resize( chemical::num_canonical_aas );
	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) P_AA_n_[ ii ].resize( 14 );

	while ( stream ) {
		using namespace ObjexxFCL::fmt;
		stream >> bite( 3, id ) >> skip( 1 ) >> bite( 2, n ) >> skip( 1 ) >> bite( 9, probability ) >> skip;
		if ( stream ) {
			assert( ( n >= 1 ) && ( n <= 14 ) ); // Support n in [1,14]
			assert( ( probability >= Probability( 0.0 ) ) && ( probability <= Probability( 1.0 ) ) );
			AA aa = aa_from_name( id );
			assert( ( aa >= 1 ) && ( aa <= num_canonical_aas ) );
			//AminoAcidKey const & key( AminoAcidKeys::key( id ) );
			P_AA_n_[ aa ][ n ] = probability;
		} //! ADD INPUT ERROR HANDLING
	}
	stream.close();

#ifndef NDEBUG
	// Check probabilities sum to ~ 1 for each (n)
	for ( int n = 1; n <= 14; ++n ) {
		Probability probability_sum( 0.0 );
		for ( Probability_AA_n::ConstIterator i = P_AA_n_.begin(), e = P_AA_n_.end(); i != e ; ++i ) {
			probability_sum += (*i)[ n ];
		}
		assert( numeric::eq_tol( probability_sum, Probability( 1.0 ), Probability( .0001 ), Probability( .0001 ) ) );
	}
#endif
}


/// @brief Read the amino acid conditional probability wrt (phi,psi) file into P_AA_pp_
///
/// @note  Only the keys present in the file are given entries
/// @note  The file entries can be in any order
/// @note  Missing entries for a present key are assigned zero
void
P_AA::read_P_AA_pp()
{
	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys::score;
	using namespace basic::options::OptionKeys::corrections::score;
	typedef  FArray2D_Probability::IR  IR; // Index range type

	// Read the probability file and load the array
	Angle phi, psi;
	std::string id;
	Probability probability;
	utility::io::izstream stream;

	// search in the local directory first
	stream.open( option[ p_aa_pp ] );
	// then database
	if ( !stream.good() ) {
		stream.close();
		basic::database::open( stream, option[ p_aa_pp ] );
	}

	if ( !stream.good() ) utility_exit_with_message( "Unable to open p_aa_pp map!" );
	P_AA_pp_.resize( chemical::num_canonical_aas );
	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
		P_AA_pp_[ ii ].dimension( IR( 0, 35 ), IR( 0, 35 ), Probability( 0.0 ) );
	}

	while ( stream ) {
		using namespace ObjexxFCL::fmt;

		stream >> bite( 4, phi ) >> skip( 1 ) >> bite( 4, psi ) >> skip( 1 )
		>> bite( 3, id ) >> skip( 17 ) >> bite( 7, probability ) >> skip;

		if ( ( stream ) ) {
			assert( ( phi >= Angle( -180.0 ) ) && ( phi <= Angle( 180.0 ) ) );
			assert( ( psi >= Angle( -180.0 ) ) && ( psi <= Angle( 180.0 ) ) );
			assert( ( probability >= Probability( 0.0 ) ) && ( probability <= Probability( 1.0 ) ) );

			AA aa = aa_from_name( id );
			assert( ( aa >= 1 ) && ( aa <= num_canonical_aas ) );

			if ( option[ p_aa_pp_nogridshift ] ) {
				int const i_phi( numeric::mod( 36 + numeric::nint(  phi / Angle( 10.0 ) ), 36 ) );
				int const i_psi( numeric::mod( 36 + numeric::nint(  psi / Angle( 10.0 ) ), 36 ) );

				if ( probability == Probability( 0.0 ) ) probability = 1e-6;
				P_AA_pp_[ aa ]( i_phi, i_psi ) = probability;
			}
			else {
				int const i_phi( numeric::mod( 36 + numeric::nint( ( phi / Angle( 10.0 ) ) - Angle( 0.5 ) ), 36 ) );
				int const i_psi( numeric::mod( 36 + numeric::nint( ( psi / Angle( 10.0 ) ) - Angle( 0.5 ) ), 36 ) );

				if ( probability == Probability( 0.0 ) ) probability = .001; //! Hack from rosetta++ except leave .001 entries alone
				P_AA_pp_[ aa ]( i_phi, i_psi ) = probability;
			}

		} //! ADD INPUT ERROR HANDLING
	}
	stream.close();

//! P_AA_pp file is NOT a proper distribution: Some (phi,psi) bins have total probabilities of zero
//! This test must be left off until the file distribution is made proper or adapted to the file semantics
//#ifndef NDEBUG
//	// Check probabilities sum to ~ 1 for each (phi,psi)
//	for ( int i_phi = 0; i_phi <= 35; ++i_phi ) {
//		for ( int i_psi = 0; i_psi <= 35; ++i_psi ) {
//			Probability probability_sum( 0.0 );
//			for ( Probability_AA_pp::ConstIterator i = P_AA_pp.begin(), e = P_AA_pp.end(); i != e ; ++i ) {
//				probability_sum += (*i)( i_phi, i_psi );
//			}
//			assert( numeric::eq_tol( probability_sum, Probability( 1.0 ), Probability( .001 ), Probability( .001 ) ) );
//		}
//	}
//#endif
}


/// @brief Probability energies from P(aa|phi,psi)
Energy
P_AA::P_AA_pp_energy( conformation::Residue const & res ) const
{
	using namespace core::chemical;
	using numeric::conversions::degrees;
	using numeric::interpolation::periodic_range::half::bilinearly_interpolated;

	AA const aa( res.aa()); //! Need to decide if/how/where to exclude NCAAs
	if ( aa > chemical::num_canonical_aas ) return 0.0;

	if ( ! res.is_terminus()  && ! res.is_virtual_residue()  )//ToDo Also exclude chainbreaks
	{ // Probabilities for this amino acid are present in files and it is not a terminus
		Angle const phi( res.mainchain_torsion( 1 ) );
		Angle const psi( res.mainchain_torsion( 2 ) );
		if ( basic::options::option[ basic::options::OptionKeys::corrections::score::p_aa_pp_nogridshift ] ) { // the format of p_aa_pp changed from using i*10+5 to i*10 as grid
			return -std::log( numeric::interpolation::periodic_range::full::bilinearly_interpolated( phi, psi, Angle( 10.0 ), 36, P_AA_pp_[ aa ] ) / P_AA_[ aa ] );
		}
		else {
			return -std::log( bilinearly_interpolated( phi, psi, Angle( 10.0 ), 36, P_AA_pp_[ aa ] ) / P_AA_[ aa ] );
		}
	} else { // Probabilities for this amino acid aren't present in files or it is a terminus
		return Energy( 0.0 );
	}
}


/// @brief Probability energies from P(aa|phi,psi): Low level calculation for non-terminus position
Energy
P_AA::P_AA_pp_energy( chemical::AA const aa, Angle const phi, Angle const psi ) const
{
	using numeric::interpolation::periodic_range::half::bilinearly_interpolated;

	if ( aa <= chemical::num_canonical_aas ) {
		// Probabilities for this amino acid are present in files
		if ( basic::options::option[ basic::options::OptionKeys::corrections::score::p_aa_pp_nogridshift ] ) { // the format of p_aa_pp changed from using i*10+5 to i*10 as grid
			return -std::log( numeric::interpolation::periodic_range::full::bilinearly_interpolated( phi, psi, Angle( 10.0 ), 36, P_AA_pp_[ aa ] ) / P_AA_[ aa ] );
		}
		else {
			return -std::log( bilinearly_interpolated( phi, psi, Angle( 10.0 ), 36, P_AA_pp_[ aa ] ) / P_AA_[ aa ] );
		}
	} else { // Probabilities for this amino acid aren't present in files or it is a terminus
		return Energy( 0.0 );
	}
}

////////////////////////////////////////////////////////////////////////////////
EnergyDerivative
P_AA::get_Paa_pp_deriv(
	conformation::Residue const & res,
	id::TorsionID const & tor_id
) const
{

	using namespace core::chemical;
	using numeric::conversions::degrees;
	using numeric::interpolation::periodic_range::half::bilinearly_interpolated;

	AA const aa( res.aa() ); //! Need to decide if/how/where to exclude NCAAs
	if ( aa > chemical::num_canonical_aas )
		return 0.0;

	/// APL ARGH!!! MAGIC NUMBERS!!!
	Size const phi_id = 1;
	Size const psi_id = 2;

	if ( ! res.is_terminus() && ( tor_id.type() == id::BB && (tor_id.torsion() == phi_id || tor_id.torsion() == psi_id )) & ! res.is_virtual_residue() ) {
		 //ToDo Also exclude chainbreaks
		// Probabilities for this amino acid are present in files and it is not a terminus
		Angle const phi( res.mainchain_torsion( phi_id ));
		Angle const psi( res.mainchain_torsion( psi_id ));
		Probability dp_dphi( 0.0 ), dp_dpsi( 0.0 );

		if ( basic::options::option[ basic::options::OptionKeys::corrections::score::p_aa_pp_nogridshift ] ) { // the format of p_aa_pp changed from using i*10+5 to i*10 as grid
			Probability const interp_p = numeric::interpolation::periodic_range::full::bilinearly_interpolated( phi, psi, Angle( 10.0 ), 36, P_AA_pp_[ aa ], dp_dphi, dp_dpsi );
			//Energy Paa_ppE = -std::log( interp_p / P_AA_[ aa ] );
			switch ( tor_id.torsion()  ) {
				case phi_id :
					return /*dlog_Paa_dphi = */ -( 1.0 / interp_p ) * dp_dphi; break;
				case psi_id :
					return /*dlog_Paa_dpsi = */ -( 1.0 / interp_p ) * dp_dpsi; break;
				default :
					return EnergyDerivative( 0.0 );
			}
		}
		else {
			Probability const interp_p = bilinearly_interpolated( phi, psi, Angle( 10.0 ), 36, P_AA_pp_[ aa ], dp_dphi, dp_dpsi );
			//Energy Paa_ppE = -std::log( interp_p / P_AA_[ aa ] );
			switch ( tor_id.torsion()  ) {
				case phi_id :
					return /*dlog_Paa_dphi = */ -( 1.0 / interp_p ) * dp_dphi; break;
				case psi_id :
					return /*dlog_Paa_dpsi = */ -( 1.0 / interp_p ) * dp_dpsi; break;
				default :
					return EnergyDerivative( 0.0 );
			}
		}
	} else { // Probabilities for this amino acid aren't present in files or it is a terminus
		return EnergyDerivative( 0.0 );
	}
}


///@brief Probability energies for P(aa)
///
///@remarks No derivative function since there are no degrees of freedom to vary for a P_AA energy like for P_AA_pp.
Energy
P_AA::P_AA_energy( conformation::Residue const & res ) const {

	using namespace core::chemical;

	AA const aa( res.aa()); //! Need to decide if/how/where to exclude NCAAs
	if ( aa > chemical::num_canonical_aas )
		return 0.0;

	return -std::log( P_AA_[ aa ] );
}


} // namespace scoring
} // namespace rosetta

