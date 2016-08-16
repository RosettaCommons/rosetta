// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ProClosureEnergy.cc
/// @brief  Proline ring closure energy method class implementation
/// @author Andrew Leaver-Fay

// Unit Headers
#include <core/scoring/methods/ProClosureEnergy.hh>
#include <core/scoring/methods/ProClosureEnergyCreator.hh>

// Package Headers
#include <core/scoring/DerivVectorPair.hh>
//#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueConnection.hh>

// Utility headers
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/dihedral_deriv.hh>

// STL Headers
#include <string>
#include <math.h>

// options
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the ProClosureEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ProClosureEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new ProClosureEnergy );
}

ScoreTypes
ProClosureEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( pro_close );
	return sts;
}


using namespace numeric::constants::d;


/// ctor
ProClosureEnergy::ProClosureEnergy() :
	parent( methods::EnergyMethodCreatorOP( new ProClosureEnergyCreator ) ),
	skip_ring_closure_(false),
	n_nv_dist_sd_( pow(basic::options::option[ basic::options::OptionKeys::score::pro_close_planar_constraint ], 2) ), // Totally fictional value.  Everywhere this is used, it's actually the square that's used.  Let's calculate the square once and only once.

	/// measured from 4745 prolines from 1.25 A and higher resolution protein structures
	trans_chi4_mean_( 176.3 * numeric::constants::d::degrees_to_radians ),
	trans_chi4_sd_( 6.0158  * numeric::constants::d::degrees_to_radians ),
	cis_chi4_mean_( -2.9105 * numeric::constants::d::degrees_to_radians ),
	cis_chi4_sd_( 5.8239    * numeric::constants::d::degrees_to_radians ),

	bbN_( "N" ),
	scNV_( "NV" ),
	scCD_( "CD" ),
	bbC_( "C" ),
	bbO_( "O")
{
	set_skip_ring_closure_from_flags();
}

ProClosureEnergy::~ProClosureEnergy()
{}

/// clone
EnergyMethodOP
ProClosureEnergy::clone() const
{
	return EnergyMethodOP( new ProClosureEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentTwoBodyEnergies
/////////////////////////////////////////////////////////////////////////////
bool
ProClosureEnergy::defines_score_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	bool res_moving_wrt_eachother
) const
{
	using namespace conformation;
	using namespace chemical;

	if ( !res_moving_wrt_eachother ) return false;

	bool const res1_is_upper( ( (rsd1.aa() == aa_pro) || (rsd1.aa() == aa_dpr) ) && rsd1.is_bonded( rsd2 ) && rsd2.has_upper_connect() && rsd2.residue_connection_partner( rsd2.upper_connect().index() ) == rsd1.seqpos() );
	bool const res2_is_upper( ( (rsd2.aa() == aa_pro) || (rsd2.aa() == aa_dpr) ) && rsd2.is_bonded( rsd1 ) && rsd1.has_upper_connect() && rsd1.residue_connection_partner( rsd1.upper_connect().index() ) == rsd2.seqpos() );

	return (res1_is_upper || res2_is_upper);
}

void
ProClosureEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	using namespace conformation;
	using namespace chemical;

	if ( rsd1.is_virtual_residue() ) return;
	if ( rsd2.is_virtual_residue() ) return;

	bool const res1_is_upper( ( (rsd1.aa() == aa_pro) || (rsd1.aa() == aa_dpr) ) && rsd1.is_bonded( rsd2 ) && rsd2.has_upper_connect() && rsd2.residue_connection_partner( rsd2.upper_connect().index() ) == rsd1.seqpos() );
	bool const res2_is_upper( ( (rsd2.aa() == aa_pro) || (rsd2.aa() == aa_dpr) ) && rsd2.is_bonded( rsd1 ) && rsd1.has_upper_connect() && rsd1.residue_connection_partner( rsd1.upper_connect().index() ) == rsd2.seqpos() );

	if ( res1_is_upper || res2_is_upper ) {

		Residue const & upper_res( res1_is_upper ? rsd1 : rsd2 );
		Residue const & lower_res( res1_is_upper ? rsd2 : rsd1 );

		Real chi4 = measure_chi4( lower_res, upper_res );
		emap[ pro_close ] += chi4E( chi4 );
	}
}

void
ProClosureEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const &,
	pose::Pose const &,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	using namespace conformation;
	using namespace chemical;

	bool const res1_is_upper( ( (rsd1.aa() == aa_pro) || (rsd1.aa() == aa_dpr) ) && rsd1.is_bonded( rsd2 ) && rsd2.has_upper_connect() && rsd2.residue_connection_partner( rsd2.upper_connect().index() ) == rsd1.seqpos() );

	conformation::Residue const & upper_res( res1_is_upper ? rsd1 : rsd2 );
	conformation::Residue const & lower_res( res1_is_upper ? rsd2 : rsd1 );

	debug_assert( (upper_res.aa() == chemical::aa_pro) || (upper_res.aa() == chemical::aa_dpr) );
	const core::Real d_multiplier = ((upper_res.aa()==chemical::aa_dpr) ? -1.0 : 1.0); //A multiplier for the derivative to invert it if this is a D-amino acid.

	utility::vector1< DerivVectorPair > & upper_res_atom_derivs( res1_is_upper ? r1_atom_derivs : r2_atom_derivs );
	utility::vector1< DerivVectorPair > & lower_res_atom_derivs( res1_is_upper ? r2_atom_derivs : r1_atom_derivs );

	/// Atoms on the upper res
	/// 1. N, CD
	/// Atoms on the lower res
	/// 2. C, O
	Size N_up_id  = upper_res.atom_index( bbN_ );
	Size CD_up_id = upper_res.atom_index( scCD_ );
	Size C_lo_id  = lower_res.atom_index( bbC_ );
	Size O_lo_id  = lower_res.atom_index( bbO_ );

	// N upper
	Vector f1( 0.0 ), f2( 0.0 );
	Real chi4( 0.0 );
	numeric::deriv::dihedral_p2_cosine_deriv(
		upper_res.xyz( CD_up_id ), upper_res.xyz( N_up_id ),
		lower_res.xyz( C_lo_id ), lower_res.xyz( O_lo_id  ),
		chi4, f1, f2 );
	chi4*=d_multiplier; //Flip for D-amino acids.
	Real deriv( weights[ pro_close ] * d_multiplier * dchi4E_dchi4( chi4 ) );

	f1 *= deriv; f2 *= deriv;
	upper_res_atom_derivs[ N_up_id ].f1() += f1;
	upper_res_atom_derivs[ N_up_id ].f2() += f2;

	/// C lower
	numeric::deriv::dihedral_p2_cosine_deriv(
		lower_res.xyz( O_lo_id  ), lower_res.xyz( C_lo_id ),
		upper_res.xyz( N_up_id ), upper_res.xyz( CD_up_id ),
		chi4, f1, f2 );
	f1 *= deriv; f2 *= deriv;
	lower_res_atom_derivs[ C_lo_id ].f1() += f1;
	lower_res_atom_derivs[ C_lo_id ].f2() += f2;


	// CD upper
	//f1 = 0.0; f2 = 0.0;
	numeric::deriv::dihedral_p1_cosine_deriv(
		upper_res.xyz( CD_up_id ), upper_res.xyz( N_up_id ),
		lower_res.xyz( C_lo_id ), lower_res.xyz( O_lo_id  ),
		chi4, f1, f2 );
	f1 *= deriv; f2 *= deriv;
	upper_res_atom_derivs[ CD_up_id ].f1() += f1;
	upper_res_atom_derivs[ CD_up_id ].f2() += f2;

	// O lower
	numeric::deriv::dihedral_p1_cosine_deriv(
		lower_res.xyz( O_lo_id  ), lower_res.xyz( C_lo_id ),
		upper_res.xyz( N_up_id ), upper_res.xyz( CD_up_id ),
		chi4, f1, f2 );
	f1 *= deriv; f2 *= deriv;
	lower_res_atom_derivs[ O_lo_id ].f1() += f1;
	lower_res_atom_derivs[ O_lo_id ].f2() += f2;
}


/// @brief Penalize the pucker-up residue type if its chi1 is positive;
/// penalize the pucker-down residue type if its chi1 is negative.  Only
/// applies this penalty when the other_residue is the next polymeric residue
/// after pro_residue (i+1), unless pro_residue is an upper_term,
/// in which case it applies the penalty for pro_residue's previous polymeric
/// residue.
/// @details Commented out.  Apparently the bump energy is not calculated for the pro_close term.
void
ProClosureEnergy::bump_energy_full(
	conformation::Residue const & ,
	conformation::Residue const & ,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const
{
}

/// @brief Penalize the pucker-up residue type if its chi1 is positive;
/// penalize the pucker-down residue type if its chi1 is negative.  Only
/// applies this penalty when the other_residue is the next polymeric residue
/// after pro_residue (i+1), unless pro_residue is an upper_term,
/// in which case it applies the penalty for pro_residue's previous polymeric
/// residue.
/// @details Commented out.  Apparently the bump energy is not calculated for the pro_close term.
void
ProClosureEnergy::bump_energy_backbone(
	conformation::Residue const & ,
	conformation::Residue const & ,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap &
) const
{
}


bool
ProClosureEnergy::defines_intrares_energy( EnergyMap const & ) const
{
	return true;
}


void
ProClosureEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap & emap
) const
{
	if ( skip_ring_closure() ) return;
	if ( (rsd.aa() == chemical::aa_pro) || (rsd.aa() == chemical::aa_dpr) ) {
		if ( rsd.is_virtual_residue() ) return;
		Distance const dist2 = rsd.xyz( bbN_ ).distance_squared( rsd.xyz( scNV_ ) );

		//Note that n_nv_dist_sd_ is the SQUARE of the standard deviation
		emap[ pro_close ] += dist2 / ( n_nv_dist_sd_ );
	}
}


bool
ProClosureEnergy::defines_intrares_energy_for_residue(
	conformation::Residue const & res
) const
{
	return ((res.aa() == chemical::aa_pro) || (res.aa() == chemical::aa_dpr));
}

void
ProClosureEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,
	pose::Pose const &,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{
	if ( skip_ring_closure() ) return;
	debug_assert ( (rsd.aa() == chemical::aa_pro) || (rsd.aa() == chemical::aa_dpr) );

	debug_assert( rsd.has( scNV_ ) );
	debug_assert( rsd.has( bbN_ ) );
	if ( rsd.is_virtual_residue() ) return;

	Size NV_ind = rsd.atom_index( scNV_ );
	Size N_ind = rsd.atom_index( bbN_ );

	Vector const & nv_pos( rsd.xyz( NV_ind ));
	Vector const & n_pos(  rsd.xyz( N_ind ));
	/// Numeric deriv version to consolidate code.

	Vector f1( 0.0 ), f2( 0.0 );
	Distance dist( 0.0 );
	numeric::deriv::distance_f1_f2_deriv( nv_pos, n_pos, dist, f1, f2 );
	Real deriv( weights[ pro_close ] * 2 * dist / ( n_nv_dist_sd_ ));
	f1 *= deriv; f2 *= deriv;

	atom_derivs[ NV_ind ].f1() += f1;
	atom_derivs[ NV_ind ].f2() += f2;
	atom_derivs[ N_ind  ].f1() -= f1;
	atom_derivs[ N_ind  ].f2() -= f2;
}

/// @brief ProClosureEnergy Energy is context independent and thus
/// indicates that no context graphs need to
/// be maintained by class Energies
void
ProClosureEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
)
const
{}

/// @brief Queries whether the user has set the -score::no_pro_close_ring_closure flag.
/// If he/she has, this sets skip_ring_closure_ to 'true'.
void ProClosureEnergy::set_skip_ring_closure_from_flags() {
	skip_ring_closure_ = basic::options::option[ basic::options::OptionKeys::score::no_pro_close_ring_closure ].user();
	return;
}

Real
ProClosureEnergy::measure_chi4(
	conformation::Residue const & lower_res,
	conformation::Residue const & upper_res
) const
{
	using namespace numeric;
	using namespace numeric::constants::d;

	debug_assert( (upper_res.aa() == chemical::aa_pro) || (upper_res.aa() == chemical::aa_dpr) );
	debug_assert( lower_res.is_bonded( upper_res ) );
	debug_assert( lower_res.has_upper_connect() && lower_res.residue_connection_partner( lower_res.upper_connect().index() ) == upper_res.seqpos() );

	Real chi4 = dihedral_radians(
		upper_res.xyz( scCD_ ),
		upper_res.xyz( bbN_ ),
		lower_res.xyz( bbC_ ),
		lower_res.xyz( bbO_ ));
	if ( upper_res.aa() == chemical::aa_dpr ) chi4*= -1.0; //invert chi4 if this is a D-pro
	if ( chi4 < -pi_over_2 ) chi4 += pi_2; // wrap
	return chi4;
}

/// @details chi4 is in radians and should be in the range between -pi_over_2 and 3/2 pi
/// If this is a D-proline, chi4 should be inverted (multiplied by -1.0) before passing to this function.
Real
ProClosureEnergy::chi4E(
	Real chi4
) const
{
	using namespace numeric::constants::d;

	if ( chi4 > pi_over_2 ) {
		Real diff = chi4 - trans_chi4_mean_;
		return diff * diff / ( trans_chi4_sd_ * trans_chi4_sd_ );
	} else {
		Real diff = chi4 - cis_chi4_mean_;
		return diff * diff / ( cis_chi4_sd_ * cis_chi4_sd_ );
	}
}


/// @details chi4 is in radians and should be in the range between -pi and pi
/// If this is a D-proline, chi4 should be inverted (multiplied by -1.0) before passing to this function,
/// and the result should also be multiplied by -1.0 before use elsewhere.
Real
ProClosureEnergy::dchi4E_dchi4(
	Real chi4
) const
{
	using namespace numeric::constants::d;
	if ( chi4 < -pi_over_2 ) chi4 += pi_2;

	if ( chi4 > pi_over_2 ) {
		Real diff = chi4 - trans_chi4_mean_;
		return 2 * diff / ( trans_chi4_sd_ * trans_chi4_sd_ );
	} else {
		Real diff = chi4 - cis_chi4_mean_;
		return 2 * diff / ( cis_chi4_sd_ * cis_chi4_sd_ );
	}
}
core::Size
ProClosureEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

