// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ProClosureEnergy.cc
/// @brief  Proline ring closure energy method class implementation
/// @author Andrew Leaver-Fay

// Unit Headers
#include <core/scoring/methods/ProClosureEnergy.hh>
#include <core/scoring/methods/ProClosureEnergyCreator.hh>

// Package Headers
#include <core/scoring/DerivVectorPair.hh>
//#include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/constraints/DihedralConstraint.hh>
// AUTO-REMOVED #include <core/scoring/func/HarmonicFunc.hh>

// Project headers
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

// Utility headers
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/dihedral_deriv.hh>

// STL Headers
#include <string>

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
	return new ProClosureEnergy;
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
	n_nv_dist_sd_( basic::options::option[ basic::options::OptionKeys::score::pro_close_planar_constraint ] ), // totally fictional

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
{}

ProClosureEnergy::~ProClosureEnergy()
{}

/// clone
EnergyMethodOP
ProClosureEnergy::clone() const
{
	return new ProClosureEnergy;
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentTwoBodyEnergies
/////////////////////////////////////////////////////////////////////////////
bool
ProClosureEnergy::defines_score_for_residue_pair(
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	bool res_moving_wrt_eachother
) const
{
	conformation::Residue const & resl( res1.seqpos() < res2.seqpos() ? res1 : res2 );
	conformation::Residue const & resu( res1.seqpos() < res2.seqpos() ? res2 : res1 );

	return res_moving_wrt_eachother && ((resu.aa() == chemical::aa_pro) || (resu.aa() == chemical::aa_dpr)) //L-pro or D-pro
		&& resl.seqpos() + 1 == resu.seqpos() && resl.is_bonded( resu );
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

	if (( ((rsd2.aa() == aa_pro) || (rsd2.aa() == aa_dpr)) && rsd1.is_bonded( rsd2 ) &&
			rsd1.seqpos() == rsd2.seqpos() - 1 ) ||
			( ((rsd1.aa() == aa_pro) || (rsd1.aa() == aa_dpr)) && rsd1.is_bonded( rsd2 ) &&
			rsd1.seqpos() - 1 == rsd2.seqpos() )) {
		Residue const & upper_res( rsd1.seqpos() > rsd2.seqpos() ? rsd1 : rsd2 );
		Residue const & lower_res( rsd1.seqpos() > rsd2.seqpos() ? rsd2 : rsd1 );

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
	bool const r1_upper( rsd1.seqpos() > rsd2.seqpos()  );
	conformation::Residue const & upper_res( r1_upper ? rsd1 : rsd2 );
	conformation::Residue const & lower_res( r1_upper ? rsd2 : rsd1 );

	assert( (upper_res.aa() == chemical::aa_pro) || (upper_res.aa() == chemical::aa_dpr) );
	const core::Real d_multiplier = ((upper_res.aa()==chemical::aa_dpr) ? -1.0 : 1.0); //A multiplier for the derivative to invert it if this is a D-amino acid.

	utility::vector1< DerivVectorPair > & upper_res_atom_derivs( r1_upper ? r1_atom_derivs : r2_atom_derivs );
	utility::vector1< DerivVectorPair > & lower_res_atom_derivs( r1_upper ? r2_atom_derivs : r1_atom_derivs );

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
void
ProClosureEnergy::bump_energy_full(
	conformation::Residue const & ,
	conformation::Residue const & ,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const
{
	/*
	static const std::string prd_name( "PRD" );
	static const std::string pru_name( "PRU" );

	if ( pro_residue.aa() == chemical::aa_pro ) {
		if ( pro_residue.is_bonded( other_residue ) &&
				( ! pro_residue.is_upper_terminus() &&
				pro_residue.seqpos() + 1 == other_residue.seqpos() )
				||
				( pro_residue.is_upper_terminus() &&
				pro_residue.seqpos() == other_residue.seqpos() + 1 ) ) {
			if ( pro_residue.name3() == prd_name && pro_residue.chi( 1 ) < 0 ) {
				emap[ pro_close ] = 10;
			} else if ( pro_residue.name3() == pru_name && pro_residue.chi( 1 ) > 0 ) {
				emap[ pro_close ] = 10;
			}
		}
	}
	*/
}

/// @brief Penalize the pucker-up residue type if its chi1 is positive;
/// penalize the pucker-down residue type if its chi1 is negative.  Only
/// applies this penalty when the other_residue is the next polymeric residue
/// after pro_residue (i+1), unless pro_residue is an upper_term,
/// in which case it applies the penalty for pro_residue's previous polymeric
/// residue.
void
ProClosureEnergy::bump_energy_backbone(
	conformation::Residue const & ,
	conformation::Residue const & ,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap &
) const
{
	//bump_energy_full( pro_residue, other_residue, pose, sfxn, emap );
}


bool
ProClosureEnergy::defines_intrares_energy( EnergyMap const & ) const
{
	return true;
}


///
void
ProClosureEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap & emap
) const
{

	if ( (rsd.aa() == chemical::aa_pro) || (rsd.aa() == chemical::aa_dpr) ) {
		if ( rsd.is_virtual_residue() )return;
		Distance const dist2 = rsd.xyz( bbN_ ).distance_squared( rsd.xyz( scNV_ ) );
		emap[ pro_close ] += dist2 / ( n_nv_dist_sd_ * n_nv_dist_sd_ );
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
	assert ( (rsd.aa() == chemical::aa_pro) || (rsd.aa() == chemical::aa_dpr) );

	//const core::Real d_multiplier = ( (rsd.aa() == chemical::aa_dpr) ? -1.0 : 1.0 ); //Multiplier for derivatives

	assert( rsd.has( scNV_ ) );
	assert( rsd.has( bbN_ ) );
	if ( rsd.is_virtual_residue() ) return;

	Size NV_ind = rsd.atom_index( scNV_ );
	Size N_ind = rsd.atom_index( bbN_ );

	Vector const & nv_pos( rsd.xyz( NV_ind ));
	Vector const & n_pos(  rsd.xyz( N_ind ));
	/// Numeric deriv version to consolidate code.
	Vector f1( 0.0 ), f2( 0.0 );
	Distance dist( 0.0 );
	numeric::deriv::distance_f1_f2_deriv( nv_pos, n_pos, dist, f1, f2 );
	Real deriv( weights[ pro_close ] * 2 * dist / ( n_nv_dist_sd_ * n_nv_dist_sd_ ));
	f1 *= deriv; f2 *= deriv;

	atom_derivs[ NV_ind ].f1() += f1;
	atom_derivs[ NV_ind ].f2() += f2;
	atom_derivs[ N_ind  ].f1() -= f1;
	atom_derivs[ N_ind  ].f2() -= f2;

}


/*void
ProClosureEnergy::eval_intrares_atom_derivative2(
	Size const atom_index,
	conformation::Residue const & rsd,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	assert ( rsd.aa() == chemical::aa_pro );

	if ( rsd.atom_index( scNV_ ) ==  atom_index ) {

		Vector const & nv_pos( rsd.xyz( atom_index ));
		Vector const & n_pos(  rsd.xyz( bbN_ ));
		/// Numeric deriv version to consolidate code.
		Vector f1( 0.0 ), f2( 0.0 );
		Distance dist( 0.0 );
		numeric::deriv::distance_f1_f2_deriv( nv_pos, n_pos, dist, f1, f2 );
		Real deriv( weights[ pro_close ] * 2 * dist / ( n_nv_dist_sd_ * n_nv_dist_sd_ ));
		F1 += deriv * f1;
		F2 += deriv * f2;

	} else if ( rsd.atom_index( bbN_ ) == atom_index ) {
		//std::cout << "evaluating N pro-closure energy derivative" << std::endl;
		Vector const & n_pos(  rsd.xyz( atom_index ));
		Vector const & nv_pos( rsd.xyz( scNV_ ));

		{ /// Scope: Distance derivative.
		Vector f1( 0.0 ), f2( 0.0 );
		Distance dist( 0.0 );
		numeric::deriv::distance_f1_f2_deriv( n_pos, nv_pos, dist, f1, f2 );
		Real deriv( weights[ pro_close ] * 2 * dist / ( n_nv_dist_sd_ * n_nv_dist_sd_ ));
		F1 += deriv * f1;
		F2 += deriv * f2;
		}

	}
}*/

/// @brief ProClosureEnergy Energy is context independent and thus
/// indicates that no context graphs need to
/// be maintained by class Energies
void
ProClosureEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
)
const
{}

Real
ProClosureEnergy::measure_chi4(
	conformation::Residue const & lower_res,
	conformation::Residue const & upper_res
) const
{
	using namespace numeric;
	using namespace numeric::constants::d;

	assert( (upper_res.aa() == chemical::aa_pro) || (upper_res.aa() == chemical::aa_dpr) );
	assert( lower_res.is_bonded( upper_res ) );
	assert( lower_res.seqpos() == upper_res.seqpos() - 1 );

	Real chi4 = dihedral_radians(
		upper_res.xyz( scCD_ ),
		upper_res.xyz( bbN_ ),
		lower_res.xyz( bbC_ ),
		lower_res.xyz( bbO_ ));
	if(upper_res.aa() == chemical::aa_dpr) chi4*= -1.0; //invert chi4 if this is a D-pro
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

