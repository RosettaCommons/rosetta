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
#include <basic/Tracer.hh>
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

#include <numeric/interpolation/spline/Bicubic_spline.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

// Utility headers
#include <utility/io/izstream.hh>

// C++ headers
#include <utility/assert.hh>

#include <utility/vector1.hh>

//MaximCode:
static thread_local basic::Tracer TR( "core.scoring.P_AA" );


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


///////////////////////////////////////////////////////////////////////////////

/// @brief Returns true if passed a core::chemical::AA corresponding to a
/// D-amino acid, and false otherwise.
bool
P_AA::is_canonical_d_aminoacid(
	core::chemical::AA const res_aa
) const {
	return core::chemical::is_canonical_D_aa(res_aa);
}

///////////////////////////////////////////////////////////////////////////////

/// @brief When passed a d-amino acid, returns the l-equivalent.  Returns
/// aa_unk otherwise.
core::chemical::AA
P_AA::get_l_equivalent(
	core::chemical::AA const d_aa
) const {
		return core::chemical::get_L_equivalent(d_aa);
}

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
	basic::database::open( stream, "scoring/score_functions/P_AA_pp/P_AA" );

	P_AA_.resize( num_canonical_aas );

	while ( stream ) {
		using namespace ObjexxFCL::format;
		stream >> bite( 3, id ) >> skip( 1 ) >> bite( 9, probability ) >> skip;
		if ( stream ) {
		debug_assert( ( probability >= Probability( 0.0 ) ) && ( probability <= Probability( 1.0 ) ) );
			AA aa = aa_from_name( id );
		debug_assert( ( aa >= 1 ) && ( aa <= num_canonical_aas ) );
			probability_sum += probability;
			P_AA_[ aa ] = probability;
		} //! ADD INPUT ERROR HANDLING
	}
	stream.close();

	// Check probabilities sum to ~ 1
debug_assert( numeric::eq_tol( probability_sum, Probability( 1.0 ), Probability( .0001 ), Probability( .0001 ) ) );
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
	basic::database::open( stream, "scoring/score_functions/P_AA_pp/P_AA_n" );

	P_AA_n_.resize( chemical::num_canonical_aas );
	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) P_AA_n_[ ii ].resize( 14 );

	while ( stream ) {
		using namespace ObjexxFCL::format;
		stream >> bite( 3, id ) >> skip( 1 ) >> bite( 2, n ) >> skip( 1 ) >> bite( 9, probability ) >> skip;
		if ( stream ) {
		debug_assert( ( n >= 1 ) && ( n <= 14 ) ); // Support n in [1,14]
		debug_assert( ( probability >= Probability( 0.0 ) ) && ( probability <= Probability( 1.0 ) ) );
			AA aa = aa_from_name( id );
		debug_assert( ( aa >= 1 ) && ( aa <= num_canonical_aas ) );
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
	debug_assert( numeric::eq_tol( probability_sum, Probability( 1.0 ), Probability( .0001 ), Probability( .0001 ) ) );
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
	using namespace basic::options::OptionKeys::corrections;
	using namespace basic::options::OptionKeys::corrections::score;
	using namespace basic::options::OptionKeys::corrections::shapovalov_lib;
	typedef  FArray2D_Probability::IR  IR; // Index range type

	// Read the probability file and load the array
	Angle phi, psi;
	std::string id;
	Probability probability;
	//MaximCode:
	Probability minusLogProbability;
	utility::io::izstream stream;

	//MaximCode
	if (option[shapovalov_lib_fixes_enable] &&
		option[shapovalov_lib::shap_p_aa_pp_enable])
	{
		std::string _smoothingRequsted = basic::options::option[ basic::options::OptionKeys::corrections::shapovalov_lib::shap_p_aa_pp_smooth_level ];
		std::string _smoothingAsKeyword = "undefined";

		if (_smoothingRequsted.compare("1")==0)
		{
			_smoothingAsKeyword = "low_smooth";
		}
		else if (_smoothingRequsted.compare("2")==0)
		{
			_smoothingAsKeyword = "high_smooth";
		}
		else
		{
			_smoothingAsKeyword = "unknown";
		}
		TR << "shapovalov_lib::shap_p_aa_pp_smooth_level of " << _smoothingRequsted << "( aka " << _smoothingAsKeyword << " )" << " got activated." << std::endl;
	}

	// search in the local directory first
	//MaximCode:
	if (!basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable]
			|| !basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_p_aa_pp_enable]) {
		stream.open( option[ p_aa_pp ] );
	}
	else {
		stream.open( option[ shap_p_aa_pp ] );
	}

	// then database
	if (!stream.good()) {
		stream.close();
		//MaximCode:
		if (!basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable]
				|| !basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_p_aa_pp_enable]) {
			basic::database::open(stream, option[p_aa_pp]);
		} else {
			basic::database::open(stream, option[shap_p_aa_pp]);
		}
	}

	if ( !stream.good() ) utility_exit_with_message( "Unable to open p_aa_pp map!" );
	P_AA_pp_.resize( chemical::num_canonical_aas );
	for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
		P_AA_pp_[ ii ].dimension( IR( 0, 35 ), IR( 0, 35 ), Probability( 0.0 ) );
	}

	while ( stream ) {
		using namespace ObjexxFCL::format;

		//MaximCode:
		if (stream.peek() == '#') {
			std::string line;
			stream.getline(line);
			continue;
		}

		//MaximCode:
		if (!basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable]
				|| !basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_p_aa_pp_enable]) {
			stream >> bite(4, phi) >> skip(1) >> bite(4, psi) >> skip(1)
					>> bite(3, id) >> skip(17) >> bite(7, probability) >> skip;
		} else {
			stream >> phi >> psi >> skip(1) >> bite(3, id) >> skip(1)
					>> probability >> minusLogProbability >> skip;
		}

		if ( ( stream ) ) {
		debug_assert( ( phi >= Angle( -180.0 ) ) && ( phi <= Angle( 180.0 ) ) );
		debug_assert( ( psi >= Angle( -180.0 ) ) && ( psi <= Angle( 180.0 ) ) );
		debug_assert( ( probability >= Probability( 0.0 ) ) && ( probability <= Probability( 1.0 ) ) );

			AA aa = aa_from_name( id );
		debug_assert( ( aa >= 1 ) && ( aa <= num_canonical_aas ) );

			//MaximCode:
			if ( option[ p_aa_pp_nogridshift ] ||
					(basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable] && basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_p_aa_pp_enable]) ) {
						int const i_phi( numeric::mod( 36 + numeric::nint( phi / Angle( 10.0 ) ), 36 ) );
						int const i_psi( numeric::mod( 36 + numeric::nint( psi / Angle( 10.0 ) ), 36 ) );

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
//		debug_assert( numeric::eq_tol( probability_sum, Probability( 1.0 ), Probability( .001 ), Probability( .001 ) ) );
//		}
//	}
//#endif


	if ( basic::options::option[ basic::options::OptionKeys::corrections::score::use_bicubic_interpolation ] ) {

		// Now prepare the bicubic spline
		using namespace numeric;
		using namespace numeric::interpolation::spline;
		P_AA_pp_energy_splines_.resize( chemical::num_canonical_aas );
		for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
			BicubicSpline paappEspline;
			MathMatrix< Real > energy_vals( 36, 36 );
			for ( Size jj = 0; jj < 36; ++jj ) {
				for ( Size kk = 0; kk < 36; ++kk ) {
					energy_vals( jj, kk ) = -std::log( P_AA_pp_[ ii ]( jj, kk ) /P_AA_[ ii ] );
				}
			}
			BorderFlag periodic_boundary[2] = { e_Periodic, e_Periodic };
			Real start_vals[2];
			if ( basic::options::option[ basic::options::OptionKeys::corrections::score::p_aa_pp_nogridshift ] ) {
				start_vals[0] = start_vals[1] = 0.0; // if this flag is on, then the PAAPP table is not shifted and aligns with the ten-degree boundaries.
			} else {
				start_vals[0] = start_vals[1] = 5.0; // otherwise, the grid is shifted by five degrees.
			}
			Real deltas[2] = {10.0, 10.0}; // grid is 10 degrees wide
			bool lincont[2] = {false,false}; //meaningless argument for a bicubic spline with periodic boundary conditions
			std::pair< Real, Real > unused[2];
			unused[0] = std::make_pair( 0.0, 0.0 );
			unused[1] = std::make_pair( 0.0, 0.0 );
			paappEspline.train( periodic_boundary, start_vals, deltas, energy_vals, lincont, unused );
			P_AA_pp_energy_splines_[ ii ] = paappEspline;
		}
	}
}


/// @brief Probability energies from P(aa|phi,psi)
/// This function handles L- and D- canonical amino acids, as well as noncanonical alpha amino acids templated on canonical alpha amino acids.
Energy
P_AA::P_AA_pp_energy( conformation::Residue const & res ) const
{
	using namespace core::chemical;
	using numeric::conversions::degrees;

	AA const aa( is_canonical_d_aminoacid(res.aa()) ? get_l_equivalent(res.aa()) : (res.backbone_aa()==aa_unk ? res.aa() : res.backbone_aa()) ); //This handles D-canonical amino acids, as well as noncanonicals templated on an L-canonical.
	if ( aa > chemical::num_canonical_aas ) return 0.0; //Excludes noncanonicals that aren't templated on a canonical.

	const core::Real d_multiplier = is_canonical_d_aminoacid(res.aa()) ? -1.0 : 1.0 ; //A multiplier that's -1 for D-amino acids and 1 for L-amino acids, used to invert phi and psi for D.

	if ( ! res.is_terminus()  && ! res.is_virtual_residue()  )//ToDo Also exclude chainbreaks
	{ // Probabilities for this amino acid are present in files and it is not a terminus
		Angle const phi( d_multiplier*res.mainchain_torsion( 1 ) );
		Angle const psi( d_multiplier*res.mainchain_torsion( 2 ) );
		//printf("P_AA_pp: res=%lu phi=%.2f psi=%.2f\n", res.seqpos(), phi, psi); fflush(stdout); //DELETE ME
		if ( basic::options::option[ basic::options::OptionKeys::corrections::score::use_bicubic_interpolation ] ) {
			return P_AA_pp_energy_splines_[ aa ].F( phi, psi );
		} else {
			if ( basic::options::option[ basic::options::OptionKeys::corrections::score::p_aa_pp_nogridshift ] ) { // the format of p_aa_pp changed from using i*10+5 to i*10 as grid
				using numeric::interpolation::periodic_range::full::bilinearly_interpolated;
				return -std::log( bilinearly_interpolated( phi, psi, Angle( 10.0 ), 36, P_AA_pp_[ aa ] ) / P_AA_[ aa ] );
			} else {
				using numeric::interpolation::periodic_range::half::bilinearly_interpolated;
				return -std::log( bilinearly_interpolated( phi, psi, Angle( 10.0 ), 36, P_AA_pp_[ aa ] ) / P_AA_[ aa ] );
			}
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

	const chemical::AA aa2 = (is_canonical_d_aminoacid(aa) ? get_l_equivalent(aa) : aa );

	if ( aa2 <= chemical::num_canonical_aas ) {
		// Probabilities for this amino acid are present in files
		if ( basic::options::option[ basic::options::OptionKeys::corrections::score::p_aa_pp_nogridshift ] ) { // the format of p_aa_pp changed from using i*10+5 to i*10 as grid
			return -std::log( numeric::interpolation::periodic_range::full::bilinearly_interpolated( phi, psi, Angle( 10.0 ), 36, P_AA_pp_[ aa2 ] ) / P_AA_[ aa2 ] );
		}
		else {
			//return -std::log( bilinearly_interpolated( phi, psi, Angle( 10.0 ), 36, P_AA_pp_[ aa ] ) / P_AA_[ aa ] );
			numeric::MathVector< Real > args(2);
			args(0) = phi;
			args(1) = psi;
			return P_AA_pp_energy_splines_[ aa2 ].F( args );
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

	AA const aa( is_canonical_d_aminoacid(res.aa()) ? get_l_equivalent(res.aa()) : (res.backbone_aa()==aa_unk ? res.aa() : res.backbone_aa()) ); //This handles D-canonical amino acids, as well as noncanonicals templated on an L-canonical.
	if ( aa > chemical::num_canonical_aas ) return 0.0; //Excludes non-templated noncanonicals.

	const core::Real d_multiplier = is_canonical_d_aminoacid(res.aa()) ? -1.0 : 1.0 ; //A multiplier that's -1 for D-amino acids and 1 for L-amino acids, used to inverte phi and psi for D.

	/// APL ARGH!!! MAGIC NUMBERS!!!
	Size const phi_id = 1;
	Size const psi_id = 2;

	if ( res.type().is_alpha_aa() && !res.is_terminus() && ( tor_id.type() == id::BB && (tor_id.torsion() == phi_id || tor_id.torsion() == psi_id )) & ! res.is_virtual_residue() ) {
		 //ToDo Also exclude chainbreaks
		// Probabilities for this amino acid are present in files and it is not a terminus
		Angle const phi( d_multiplier*res.mainchain_torsion( phi_id ));
		Angle const psi( d_multiplier*res.mainchain_torsion( psi_id ));
		Probability dp_dphi( 0.0 ), dp_dpsi( 0.0 );

		if ( basic::options::option[ basic::options::OptionKeys::corrections::score::use_bicubic_interpolation ] ) {
			switch ( tor_id.torsion()  ) {
			case phi_id :
				return d_multiplier*P_AA_pp_energy_splines_[ aa ].dFdx( phi, psi );
			case psi_id :
				return d_multiplier*P_AA_pp_energy_splines_[ aa ].dFdy( phi, psi );
			default :
				return EnergyDerivative( 0.0 );
			}
		} else {
			if ( basic::options::option[ basic::options::OptionKeys::corrections::score::p_aa_pp_nogridshift ] ) { // the format of p_aa_pp changed from using i*10+5 to i*10 as grid
				Probability const interp_p = numeric::interpolation::periodic_range::full::bilinearly_interpolated( phi, psi, Angle( 10.0 ), 36, P_AA_pp_[ aa ], dp_dphi, dp_dpsi );
				//Energy Paa_ppE = -std::log( interp_p / P_AA_[ aa ] );
				switch ( tor_id.torsion()  ) {
				case phi_id :
					return /*dlog_Paa_dphi = */ -( 1.0 / interp_p ) * d_multiplier * dp_dphi; break;
				case psi_id :
					return /*dlog_Paa_dpsi = */ -( 1.0 / interp_p ) * d_multiplier * dp_dpsi; break;
				default :
					return EnergyDerivative( 0.0 );
				}
			} else {
				Real const interp_p = bilinearly_interpolated( phi, psi, Angle( 10.0 ), 36, P_AA_pp_[ aa ], dp_dphi, dp_dpsi );
				switch ( tor_id.torsion()  ) {
					case phi_id :
						return /*dlog_Paa_dphi = */ -( 1.0 / interp_p ) * d_multiplier * dp_dphi; break;
					case psi_id :
						return /*dlog_Paa_dpsi = */ -( 1.0 / interp_p ) * d_multiplier * dp_dpsi; break;
					default :
						return EnergyDerivative( 0.0 );
				}
			}
		}
	} else { // Probabilities for this amino acid aren't present in files or it is a terminus
		return EnergyDerivative( 0.0 );
	}
}


/// @brief Probability energies for P(aa)
///
/// @remarks No derivative function since there are no degrees of freedom to vary for a P_AA energy like for P_AA_pp.
Energy
P_AA::P_AA_energy( conformation::Residue const & res ) const {

	using namespace core::chemical;

	AA const aa( is_canonical_d_aminoacid(res.aa()) ? get_l_equivalent(res.aa()) : (res.backbone_aa()==aa_unk ? res.aa() : res.backbone_aa()) ); //This handles D-canonical amino acids, as well as noncanonicals templated on an L-canonical.
	if ( aa > chemical::num_canonical_aas )
		return 0.0;

	return -std::log( P_AA_[ aa ] );
}


} // namespace scoring
} // namespace rosetta

