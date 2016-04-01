// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Ramachandran2B.cc
/// @brief  Neighbor-dependent Ramachandran potential class implementation
/// @author Guoli Wang
/// @author Amelie Stein (amelie.stein@ucsf.edu) Oct 2012 -- rama2b lookup table for loop modeling
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/Ramachandran2B.hh>

// Package Headers
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/interpolation/periodic_range/half/interpolation.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>

// Utility Headers
#include <utility/io/izstream.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

// AS -- to get access to get_torsion_bin()
#include <core/conformation/ppo_torsion_bin.hh>

// option key includes
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2A.hh>


using namespace ObjexxFCL;

namespace core {
namespace scoring {


// @brief Auto-generated virtual destructor
Ramachandran2B::~Ramachandran2B() {}

static THREAD_LOCAL basic::Tracer T( "core.scoring.Ramachandran2B" );

Real const Ramachandran2B::binw_( 10.0 );

// AS
// only sample torsions with Rama prob above this value -- note that values are
// directly copied from the Ramachandran.cc implementation, might need tweaking
Real const Ramachandran2B::rama_sampling_thold_(0.00075 );

Ramachandran2B::Ramachandran2B() :
	ram_energ_( n_phi_, n_psi_, n_aa_, 0.0 ),
	ram_entropy_( n_aa_, 0.0 ),
	ram_energ_left_( n_phi_, n_psi_, n_aa_, n_aa_, 0.0 ),
	ram_entropy_left_( n_aa_, n_aa_, 0.0 ),
	ram_energ_right_( n_phi_, n_psi_, n_aa_, n_aa_, 0.0 ),
	ram_entropy_right_( n_aa_, n_aa_, 0.0 ),
	rama_score_limit_( 20 ),
	left_ram_probabil_( n_phi_, n_psi_, n_aa_, n_aa_, 0.0 ),
	right_ram_probabil_( n_phi_, n_psi_, n_aa_, n_aa_, 0.0 ),
	left_cdf_( n_aa_, n_aa_ ),
	left_cdf_by_torsion_bin_( n_aa_, n_aa_, conformation::n_ppo_torsion_bins ),
	n_valid_left_pp_bins_by_ppo_torbin_( n_aa_, n_aa_, conformation::n_ppo_torsion_bins ),
	right_cdf_( n_aa_, n_aa_ ),
	right_cdf_by_torsion_bin_( n_aa_, n_aa_, conformation::n_ppo_torsion_bins ),
	n_valid_right_pp_bins_by_ppo_torbin_( n_aa_, n_aa_, conformation::n_ppo_torsion_bins )
{
	read_rama();
}


///////////////////////////////////////////////////////////////////////////////

/// @brief evaluate rama score for each (protein) residue and store that score
/// in the pose.energies() object
void
Ramachandran2B::eval_rama_score_all(
	pose::Pose & pose,
	ScoreFunction const & scorefxn
) const {
	if ( scorefxn.has_zero_weight( rama ) ) return; // unnecessary, righ?

	// in pose mode, we use fold_tree.cutpoint info to exclude terminus
	// residues from rama calculation. A cutpoint could be either an artificial
	// cutpoint such as loop cutpoint or a real physical chain break such as
	// multiple-chain complex. For the artificial cutpoint, we may need to
	// calculate rama scores for cutpoint residues, but for the real chain break
	// cutpoint, we don't want to do that. So here we first loop over all the
	// residue in the protein and exclude those ones which are the cutpoints.
	// Then we loop over the cutpoint residues and add rama score for residues
	// at artificial cutpoints, i.e., cut_weight != 0.0, which means that
	// jmp_chainbreak_score is also calculated for this cutpoint. Note that the
	// default value for cut_weight here is dependent on whether
	// jmp_chainbreak_weight is set. This is to ensure that rama score for
	// termini residues are not calculated when jmp_chainbreak_weight is 0.0,
	// e.g normal pose docking.

	int const total_residue = pose.total_residue();

	// exclude chain breaks
	Energies & pose_energies( pose.energies() );

	for ( int ii = 1; ii <= total_residue; ++ii ) {
		if ( !pose.residue(ii).is_protein() || pose.residue(ii).is_terminus() ) continue;
		
		Real rama_score,dphi,dpsi;
		eval_rama_score_residue(pose.residue(ii), pose.residue(ii-1).aa(), pose.residue(ii+1).aa(), rama_score, dphi, dpsi);
		T << "Rama:eval_all: residue " << ii << " " << pose.residue(ii).name() <<
		" " << ii-1 << " " << pose.residue(ii-1).name() << " " << ii+1 << " " <<
		pose.residue(ii+1).name() << " = " << rama_score << std::endl;
		pose_energies.onebody_energies( ii )[rama] = rama_score;
	}
}


///////////////////////////////////////////////////////////////////////////////
void
Ramachandran2B::write_rama_score_all( Pose const & /*pose*/ ) const
{}


///////////////////////////////////////////////////////////////////////////////
void
Ramachandran2B::eval_rama_score_residue(
	conformation::Residue const & rsd,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const {
	using namespace numeric;

	debug_assert( rsd.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2)));

	// amw replacing anything that would set rama to 0 because of an incidental phi/psi of 0
	// (rare but possible)
	if ( rsd.type().has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT )
			|| rsd.type().has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) || rsd.is_terminus() ) { // begin or end of chain
		rama = 0.0;
		drama_dphi = 0.0;
		drama_dpsi = 0.0;
		return;
	}

	eval_rama_score_residue( rsd.aa(), phi, psi, rama, drama_dphi, drama_dpsi );
}

///////////////////////////////////////////////////////////////////////////////
// modified by GL
void
Ramachandran2B::eval_rama_score_residue(
	conformation::Residue const &center,
	chemical::AA const left_aa,
	chemical::AA const right_aa,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const {
	using namespace numeric;

	debug_assert( center.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( center.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( center.mainchain_torsion(2)));

	if ( center.type().has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT )
			|| center.type().has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) ) { //phi == 0.0 || psi == 0.0 ) { // begin or end of chain
		rama = 0.0;
		drama_dphi = 0.0;
		drama_dpsi = 0.0;
		return;
	}

	if ( ! basic::options::option[ basic::options::OptionKeys::score::ramaneighbors ] ) {
		eval_rama_score_residue( center.aa(), phi, psi, rama, drama_dphi, drama_dpsi );
	} else {
		rama = eval_rama_score_residue( phi, psi, center.aa(), left_aa, right_aa, drama_dphi, drama_dpsi );
	}
}

Real
Ramachandran2B::eval_rama_score_residue(
	Real phi,
	Real psi,
	chemical::AA const res,
	chemical::AA const left_aa,
	chemical::AA const right_aa
) const {
	Real drama_dphi, drama_dpsi;
	return eval_rama_score_residue( phi, psi, res, left_aa, right_aa, drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::eval_rama_score_residue(
	Real phi,
	Real psi,
	chemical::AA const res,
	chemical::AA const left_aa,
	chemical::AA const right_aa,
	Real & drama_dphi,
	Real & drama_dpsi
) const {
	Real rama_L( 0.0 ), drama_dphi_L( 0.0 ), drama_dpsi_L( 0.0 );
	Real rama_R( 0.0 ), drama_dphi_R( 0.0 ), drama_dpsi_R( 0.0 );
	Real rama_0( 0.0 ), drama_dphi_0( 0.0 ), drama_dpsi_0( 0.0 );
	rama_L = RamaE_Lower( phi, psi, res, left_aa, drama_dphi_L, drama_dpsi_L );
	rama_R = RamaE_Upper( phi, psi, res, right_aa, drama_dphi_R, drama_dpsi_R );
	rama_0 = RamaE( phi, psi, res, drama_dphi_0, drama_dpsi_0 );

	drama_dphi = drama_dphi_L + drama_dphi_R - drama_dphi_0;
	drama_dpsi = drama_dpsi_L + drama_dpsi_R - drama_dpsi_0;

	return rama_L + rama_R - rama_0; // rama
}

void
Ramachandran2B::IdealizeRamaEnergy(
	Real const phi,
	Real const psi,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi,
	Real const entropy,
	FArray2A< Real > const & rama_for_res
) const {
	using namespace numeric::interpolation::periodic_range::half;
	Real interp_E = bilinearly_interpolated( phi, psi, binw_, n_phi_, rama_for_res, drama_dphi, drama_dpsi );
	rama = entropy + interp_E;
	// std::cout << "Rama::eval_res: " <<  interp_E << " rama " << rama << std::endl;

	if ( ! basic::options::option[basic::options::OptionKeys::corrections::score::rama_not_squared] ) {
		if ( rama > 1.0 ) {
			Real rama_squared = rama * rama;
			if ( rama_squared > rama_score_limit_ ) {
				drama_dphi = 0.0;
				drama_dpsi = 0.0;
				rama = rama_score_limit_;
			} else {
				drama_dphi *= 2.0 * rama;
				drama_dpsi *= 2.0 * rama;
				rama = rama_squared;
			}
		}
	}
	// std::cout << " rama: " << rama << " dphi " << drama_dphi << " dpsi " << drama_dpsi << std::endl;
}

// end modification

///////////////////////////////////////////////////////////////////////////////
// modified by GL according to Andrew's suggestion
Real
Ramachandran2B::RamaE_Lower(
	conformation::Residue const &rsd,
	chemical::AA const &neighbor
) const {
	Real drama_dphi, drama_dpsi;
	return RamaE_Lower( rsd, neighbor, drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::RamaE_Lower(
	conformation::Residue const &rsd,
	chemical::AA const &neighbor,
	Real & drama_dphi,
	Real & drama_dpsi
) const {
	using namespace numeric;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	debug_assert( rsd.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2)));

	if ( rsd.type().has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT )
			|| rsd.type().has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) ) { //phi == 0.0 || psi == 0.0 ) { // begin or end of chain
		return 0.0;
	}

	// if neighbor independent protocol is selected, return 0.0
	if ( ! option[ score::ramaneighbors ] ) {
		return 0.0;
	}

	return RamaE_Lower( phi, psi, rsd.aa(), neighbor, drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::RamaE_Lower(
	Real phi,
	Real psi,
	chemical::AA const & rsd,
	chemical::AA const & neighbor
) const {
	Real drama_dphi, drama_dpsi;
	return RamaE_Lower( phi, psi, rsd, neighbor, drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::RamaE_Lower(
	Real phi,
	Real psi,
	chemical::AA const & rsd,
	chemical::AA const & neighbor,
	Real & drama_dphi,
	Real & drama_dpsi
) const {
	FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
	FArray2A< Real > const rama_for_res( ram_energ_left_(1, 1, rsd, neighbor), zero_index, zero_index );
	Real entropy = ram_entropy_left_(rsd, neighbor);

	Real rama;
	IdealizeRamaEnergy( phi, psi, rama, drama_dphi, drama_dpsi, entropy, rama_for_res );
	return rama;
}

Real
Ramachandran2B::RamaE_Upper(
	conformation::Residue const &rsd,
	chemical::AA const &neighbor
) const {
	Real drama_dphi, drama_dpsi;
	return RamaE_Upper( rsd, neighbor, drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::RamaE_Upper(
	conformation::Residue const &rsd,
	chemical::AA const &neighbor,
	Real & drama_dphi,
	Real & drama_dpsi
) const {
	using namespace numeric;

	debug_assert( rsd.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2)));

	if ( rsd.type().has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT )
			|| rsd.type().has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) ) { // begin or end of chain
		return 0.0;
	}

	// if neighbor independent protocol is selected, return 0.0
	if ( ! basic::options::option[ basic::options::OptionKeys::score::ramaneighbors ] ) {
		return 0.0;
	}

	return RamaE_Upper( phi, psi, rsd.aa(), neighbor, drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::RamaE_Upper(
	Real phi,
	Real psi,
	chemical::AA const & rsd,
	chemical::AA const & neighbor
) const {
	Real drama_dphi, drama_dpsi;
	return RamaE_Upper( phi, psi, rsd, neighbor, drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::RamaE_Upper(
	Real phi,
	Real psi,
	chemical::AA const & rsd,
	chemical::AA const & neighbor,
	Real & drama_dphi,
	Real & drama_dpsi
) const {
	FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1 );
	FArray2A< Real > const rama_for_res( ram_energ_right_( 1, 1, rsd, neighbor ), zero_index, zero_index );
	Real entropy = ram_entropy_right_( rsd, neighbor );

	Real rama;
	IdealizeRamaEnergy( phi, psi, rama, drama_dphi, drama_dpsi, entropy, rama_for_res );
	return rama;
}

Real
Ramachandran2B::RamaE(
	conformation::Residue const &rsd
) const {
	Real drama_dphi(0.0), drama_dpsi(0.0);
	return RamaE( rsd, drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::RamaE(
	conformation::Residue const &rsd,
	Real & drama_dphi,
	Real & drama_dpsi
) const {
	using namespace numeric;

	debug_assert( rsd.is_protein() );

	Real const phi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1)));
	Real const psi
		( nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2)));

	if ( rsd.type().has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT )
			|| rsd.type().has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) ) { //phi == 0.0 || psi == 0.0 ) { // begin or end of chain
		return 0.0;
	}

	return RamaE( phi, psi, rsd.aa(), drama_dphi, drama_dpsi );
}

Real
Ramachandran2B::RamaE(
	Real phi,
	Real psi,
	chemical::AA const & rsd,
	Real &drama_dphi,
	Real &drama_dpsi
) const {
	Real ramaE(0.0);
	eval_rama_score_residue( rsd, phi, psi, ramaE, drama_dphi, drama_dpsi );
	return ramaE;
}

///////////////////////////////////////////////////////////////////////////////
///
Real
Ramachandran2B::eval_rama_score_residue(
	AA const res_aa,
	Real const phi,
	Real const psi
) const {

	Real rama, drama_dphi, drama_dpsi;
	eval_rama_score_residue( res_aa, phi, psi, rama, drama_dphi, drama_dpsi );
	return rama;
}

///////////////////////////////////////////////////////////////////////////////
///
void
Ramachandran2B::eval_rama_score_residue(
	AA const res_aa,
	Real const phi,
	Real const psi,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const {
	using namespace numeric;

	FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
	FArray2A< Real > const rama_for_res( ram_energ_(1, 1, res_aa ), zero_index, zero_index );
	Real entropy = ram_entropy_(res_aa);

	IdealizeRamaEnergy( phi, psi, rama, drama_dphi, drama_dpsi, entropy, rama_for_res );
}


///////////////////////////////////////////////////////////////////////////////
/// Sample phi/psi torsions with probabilities proportionate to their
/// Ramachandran probabilities
/// Note -- this function had previously required that the option
/// loops::nonpivot_torsion_sampling be active.  This function now
/// performs a just-in-time check to initialize these tables the first
/// time they are requested -- To properly multi-thread this code, the
/// function should nab a mutex so that no two threads try to execute
/// the code at once.
void
Ramachandran2B::random_phipsi_from_rama_left(
	AA const left_aa,
	AA const pos_aa,
	Real & phi,
	Real & psi
) const {
	draw_random_phi_psi_from_cdf( left_cdf_( left_aa, pos_aa ), phi, psi );
}


///////////////////////////////////////////////////////////////////////////////
/// @details Sample phi/psi torsions with probabilities proportionate to their
/// Ramachandran probabilities
/// Note -- this function had previously required that the option
/// loops::nonpivot_torsion_sampling be active.  This function now
/// performs a just-in-time check to initialize these tables the first
/// time they are requested -- To properly multi-thread this code, the
/// function should nab a mutex so that no two threads try to execute
/// the code at once.
void
Ramachandran2B::random_phipsi_from_rama_right(
	AA const pos_aa,
	AA const right_aa,
	Real & phi,
	Real & psi
) const {
	draw_random_phi_psi_from_cdf( right_cdf_( right_aa, pos_aa ), phi, psi );
}


///////////////////////////////////////////////////////////////////////////////
/// Sample phi/psi torsions with probabilities proportionate to their
/// Ramachandran probabilities -- this version performs lookup restricted to specified torsion bins
/// based on random_phipsi_from_rama and has the same issue for parallel running

/// @author Amelie Stein (amelie.stein@ucsf.edu)
/// @date Fri May 11 15:52:01 PDT 2012
/// @details returns a random phi/psi combination within the given torsion bin -- WARNING: this will only work for the torsion bins that are currently implemented
void
Ramachandran2B::random_phipsi_from_rama_by_torsion_bin_left(
	AA const left_aa,
	AA const pos_aa,
	Real & phi,
	Real & psi,
	conformation::ppo_torsion_bin const torsion_bin
) const {
	draw_random_phi_psi_from_cdf( left_cdf_by_torsion_bin_( left_aa, pos_aa, torsion_bin ), phi, psi );
} // random_phipsi_from_rama_by_torsion_bin_left

void
Ramachandran2B::random_phipsi_from_rama_by_torsion_bin_right(
	AA const pos_aa,
	AA const right_aa,
	Real & phi,
	Real & psi,
	conformation::ppo_torsion_bin const torsion_bin
) const {
	draw_random_phi_psi_from_cdf( right_cdf_by_torsion_bin_( right_aa, pos_aa, torsion_bin ), phi, psi );
}

void
Ramachandran2B::get_entries_per_torsion_bin_left(
	AA const left_aa,
	AA const pos_aa,
	std::map< conformation::ppo_torsion_bin, core::Size > & tb_frequencies
) const {
	tb_frequencies[conformation::ppo_torbin_A] = n_valid_left_pp_bins_by_ppo_torbin_( left_aa, pos_aa, conformation::ppo_torbin_A );
	tb_frequencies[conformation::ppo_torbin_B] = n_valid_left_pp_bins_by_ppo_torbin_( left_aa, pos_aa, conformation::ppo_torbin_B );
	tb_frequencies[conformation::ppo_torbin_E] = n_valid_left_pp_bins_by_ppo_torbin_( left_aa, pos_aa, conformation::ppo_torbin_E );
	tb_frequencies[conformation::ppo_torbin_G] = n_valid_left_pp_bins_by_ppo_torbin_( left_aa, pos_aa, conformation::ppo_torbin_G );
	tb_frequencies[conformation::ppo_torbin_X] = n_valid_left_pp_bins_by_ppo_torbin_( left_aa, pos_aa, conformation::ppo_torbin_X );
}


void
Ramachandran2B::get_entries_per_torsion_bin_right(
	AA const pos_aa,
	AA const right_aa,
	std::map< conformation::ppo_torsion_bin, core::Size > & tb_frequencies
) const {
	tb_frequencies[conformation::ppo_torbin_A] = n_valid_right_pp_bins_by_ppo_torbin_( right_aa, pos_aa, conformation::ppo_torbin_A );
	tb_frequencies[conformation::ppo_torbin_B] = n_valid_right_pp_bins_by_ppo_torbin_( right_aa, pos_aa, conformation::ppo_torbin_B );
	tb_frequencies[conformation::ppo_torbin_E] = n_valid_right_pp_bins_by_ppo_torbin_( right_aa, pos_aa, conformation::ppo_torbin_E );
	tb_frequencies[conformation::ppo_torbin_G] = n_valid_right_pp_bins_by_ppo_torbin_( right_aa, pos_aa, conformation::ppo_torbin_G );
	tb_frequencies[conformation::ppo_torbin_X] = n_valid_right_pp_bins_by_ppo_torbin_( right_aa, pos_aa, conformation::ppo_torbin_X );
}

Size Ramachandran2B::n_phi_bins() const { return n_phi_; }

Size Ramachandran2B::n_psi_bins() const { return n_psi_; }

Real
Ramachandran2B::rama_bin_probability_left(
	core::chemical::AA aa_left,
	core::chemical::AA aa_center,
	Real phi,
	Real psi
) const {
	phi = numeric::nonnegative_principal_angle_degrees( phi );
	psi = numeric::nonnegative_principal_angle_degrees( psi );
	Size phi_ind = static_cast< Size > ( floor( phi / binw_ ) ) + 1;
	Size psi_ind = static_cast< Size > ( floor( psi / binw_ ) ) + 1;
	return left_ram_probabil_( phi_ind, psi_ind, aa_left, aa_center );
}

Real
Ramachandran2B::rama_bin_probability_right(
	core::chemical::AA aa_center,
	core::chemical::AA aa_right,
	Real phi,
	Real psi
) const {
	phi = numeric::nonnegative_principal_angle_degrees( phi );
	psi = numeric::nonnegative_principal_angle_degrees( psi );
	Size phi_ind = static_cast< Size > ( floor( phi / binw_ ) ) + 1;
	Size psi_ind = static_cast< Size > ( floor( psi / binw_ ) ) + 1;
	return right_ram_probabil_( phi_ind, psi_ind, aa_right, aa_center );
}

Real Ramachandran2B::minimum_sampling_probability() const {
	return rama_sampling_thold_;
}

utility::vector1< Real > const &
Ramachandran2B::left_cdf(
	core::chemical::AA aa_left,
	core::chemical::AA aa_center
) const {
	return left_cdf_( aa_left, aa_center );
}

utility::vector1< Real > const &
Ramachandran2B::right_cdf(
	core::chemical::AA aa_center,
	core::chemical::AA aa_right
) const {
	return right_cdf_( aa_right, aa_center );
}

utility::vector1< Real > const &
Ramachandran2B::left_cdf_for_torsion_bin(
	chemical::AA aa_left,
	chemical::AA aa_center,
	conformation::ppo_torsion_bin torsion_bin
) const {
	return left_cdf_by_torsion_bin_( aa_left, aa_center, torsion_bin );
}

utility::vector1< Real > const &
Ramachandran2B::right_cdf_for_torsion_bin(
	chemical::AA aa_center,
	chemical::AA aa_right,
	conformation::ppo_torsion_bin torsion_bin
) const {
	return right_cdf_by_torsion_bin_( aa_right, aa_center, torsion_bin );
}

// Guoli Wang
void
Ramachandran2B::read_rama()
{
	using namespace basic::options;

	using namespace basic::options::OptionKeys;

	int aa_num( 0 ), aa_num_left( 0 ), aa_num_right( 0 );
	int phi_bin( 0 ), psi_bin( 0 ), ss_type( 0 );
	int tCounts( 0 );
	Real tProb( 0.0 ), tEnergy( 0.0 );

	Size line_count( 0 );

	//utility::io::izstream  iunit;
#ifndef WIN32
#ifndef __CYGWIN__
	clock_t starttime = clock();
#endif
#endif
	std::string energyFileName = basic::options::option[ in::file::rama2b_map ]().name() ; // "wrapHDPprobs36.both";
	T << "Read in ramachandran map: " <<  energyFileName << std::endl;
	utility::io::izstream iRamaEnergy;
	basic::database::open( iRamaEnergy, energyFileName );
	while ( ! iRamaEnergy.eof() ) {
		++line_count;
		iRamaEnergy >> aa_num >> aa_num_left >> aa_num_right >> ss_type >> phi_bin >> psi_bin >> tCounts >> tProb >> tEnergy;
		// std::cout << " aa_num " << aa_num << " aa_num_left " << aa_num_left << " aa_num_right " << aa_num_right << " ss_type " << ss_type <<
		//    " phi_bin " << phi_bin << " psi_bin " << psi_bin << " tProb " << tProb << " tEnergy " << tEnergy << std::endl;
		if ( aa_num > n_aa_ ) continue;

		int phiIndex = phi_bin / 10 + 1;
		int psiIndex = psi_bin / 10 + 1;
		Real entropy = -1.0 * tProb * tEnergy;

		if ( aa_num_left == nullaa && aa_num_right == nullaa ) {
			ram_energ_( phiIndex, psiIndex, aa_num ) = tEnergy;
			ram_entropy_( aa_num ) += entropy;
		} else if ( aa_num_left != nullaa ) {
			ram_energ_left_( phiIndex, psiIndex, aa_num, aa_num_left ) = tEnergy;
			ram_entropy_left_( aa_num, aa_num_left ) += entropy;
			left_ram_probabil_( phiIndex, psiIndex, aa_num_left, aa_num ) = tProb;
		} else if ( aa_num_right != nullaa ) {
			ram_energ_right_( phiIndex, psiIndex, aa_num, aa_num_right ) = tEnergy;
			ram_entropy_right_( aa_num, aa_num_right ) += entropy;
			right_ram_probabil_( phiIndex, psiIndex, aa_num_right, aa_num ) = tProb; // APL changing the indexing here to be consistent w/ left_ram_probabil_
		}
	}

	iRamaEnergy.close();

#ifndef WIN32
#ifndef __CYGWIN__
	clock_t stoptime = clock();
	T << "Reading Rama from database took " << ((double) stoptime - starttime)/CLOCKS_PER_SEC << " seconds" << std::endl;
#endif
#endif

	initialize_rama_sampling_tables();
}

void
Ramachandran2B::initialize_rama_sampling_tables()
{
	init_rama_sampling_table( conformation::ppo_torbin_X, left_ram_probabil_, left_cdf_, left_cdf_by_torsion_bin_, n_valid_left_pp_bins_by_ppo_torbin_ );
	init_rama_sampling_table( conformation::ppo_torbin_A, left_ram_probabil_, left_cdf_, left_cdf_by_torsion_bin_, n_valid_left_pp_bins_by_ppo_torbin_ );
	init_rama_sampling_table( conformation::ppo_torbin_B, left_ram_probabil_, left_cdf_, left_cdf_by_torsion_bin_, n_valid_left_pp_bins_by_ppo_torbin_ );
	init_rama_sampling_table( conformation::ppo_torbin_E, left_ram_probabil_, left_cdf_, left_cdf_by_torsion_bin_, n_valid_left_pp_bins_by_ppo_torbin_ );
	init_rama_sampling_table( conformation::ppo_torbin_G, left_ram_probabil_, left_cdf_, left_cdf_by_torsion_bin_, n_valid_left_pp_bins_by_ppo_torbin_ );

	init_rama_sampling_table( conformation::ppo_torbin_X, right_ram_probabil_, right_cdf_, right_cdf_by_torsion_bin_, n_valid_right_pp_bins_by_ppo_torbin_ );
	init_rama_sampling_table( conformation::ppo_torbin_A, right_ram_probabil_, right_cdf_, right_cdf_by_torsion_bin_, n_valid_right_pp_bins_by_ppo_torbin_ );
	init_rama_sampling_table( conformation::ppo_torbin_B, right_ram_probabil_, right_cdf_, right_cdf_by_torsion_bin_, n_valid_right_pp_bins_by_ppo_torbin_ );
	init_rama_sampling_table( conformation::ppo_torbin_E, right_ram_probabil_, right_cdf_, right_cdf_by_torsion_bin_, n_valid_right_pp_bins_by_ppo_torbin_ );
	init_rama_sampling_table( conformation::ppo_torbin_G, right_ram_probabil_, right_cdf_, right_cdf_by_torsion_bin_, n_valid_right_pp_bins_by_ppo_torbin_ );

}


/// @details Compute the cumulative distribution function for a particular
/// phi/psi/omega torsion bin over all amino-acid pairs and store the result in
/// the cdf (so long as torsion_bin is ppo_torbin_X) and the cdf_by_torsion_bin
/// talbes, and store a count of the number of valid phi/psi bins in the
/// n_valid_pp_bins_by_ppo_torbin table.  A phi/psi bin is valid if its
/// probability is greater than rama_sampling_thold_ and if it falls within
/// the input torsion_bin (as answered by core::conformation::get_torsion_bin)
void
Ramachandran2B::init_rama_sampling_table(
	const conformation::ppo_torsion_bin torsion_bin,
	ObjexxFCL::FArray4D< Real > const & ram_probability,
	ObjexxFCL::FArray2D< utility::vector1< Real > > & cdf,
	ObjexxFCL::FArray3D< utility::vector1< Real > > & cdf_by_torsion_bin,
	ObjexxFCL::FArray3D< Size > & n_valid_pp_bins_by_ppo_torbin
) const {
	FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
	utility::vector1< Real > inner_cdf( n_phi_ * n_psi_, 0.0 );
	for ( int ii=1; ii <= n_aa_; ++ii ) {
		for ( int jj=1; jj <= n_aa_; ++jj ) {

			std::fill( inner_cdf.begin(), inner_cdf.end(), Real( 0.0 ) );

			FArray2A< Real > const iijj_rama_prob( ram_probability(1, 1, jj, ii), zero_index, zero_index );
			Size actual_allowed = 0;

			// current_rama_sampling_table[left_aa][aa][right_aa].resize(max_allowed); // I think this is resized later anyway
			Real allowed_probability_sum = 0.0;
			for ( int kk = 0; kk < (int) n_phi_; ++kk ) {
				for ( int ll = 0; ll < (int) n_psi_; ++ll ) {

					// store the cumulative sum up to the current table entry
					inner_cdf[ kk*n_psi_ + ll +  1 ] = allowed_probability_sum;

					Real kkll_prob = iijj_rama_prob(kk,ll);
					if ( kkll_prob < rama_sampling_thold_ ) continue;

					// from kk and ll, compute the phi and psi values within the range [-180,180]
					Real const cur_phi = binw_ * ( kk - ( kk > n_phi_ / 2 ? n_phi_ : 0 ));
					Real const cur_psi = binw_ * ( ll - ( ll > n_psi_ / 2 ? n_psi_ : 0 ));

					conformation::ppo_torsion_bin cur_tb = conformation::ppo_torbin_X;
					if ( torsion_bin != conformation::ppo_torbin_X ) {
						//  AS -- how can we get the factor properly / without hard-coding? - also: this takes very long...
						cur_tb = core::conformation::get_torsion_bin(cur_phi, cur_psi);
					}

					if ( cur_tb == torsion_bin ) {
						++actual_allowed;
						allowed_probability_sum += kkll_prob;
					}
				}
			}

			// now normalize the cdf vector by scaling by 1/allowed_probability_sum so that it represents a CDF that sums to 1.
			// and mark all zero-probability bins with the CDF value of their predecessor.
			Real const inv_allowed_probability_sum = 1 / allowed_probability_sum;
			for ( Size kk = 1; kk <= n_phi_ * n_psi_; ++kk ) {
				inner_cdf[ kk ] *= inv_allowed_probability_sum;
				if ( inner_cdf[ kk ] == 0 && kk != 1 ) {
					inner_cdf[ kk ] = inner_cdf[ kk-1 ];
				}
			}

			// now update the cdf, cdf_by_torsion_bin, and n_valid_pp_bins_by_ppo_torbin tables
			if ( torsion_bin == conformation::ppo_torbin_X ) {
				cdf( jj, ii ) = inner_cdf;
			}
			cdf_by_torsion_bin( jj, ii, torsion_bin ) = inner_cdf;
			n_valid_pp_bins_by_ppo_torbin( jj, ii, torsion_bin ) = actual_allowed;
		}
	}
}

void
Ramachandran2B::draw_random_phi_psi_from_cdf(
	utility::vector1< Real > const & cdf,
	Real & phi,
	Real & psi
) const {
	// the bin index can be unpacked to give the phi and psi indices
	Size bin_from_cdf = numeric::random::pick_random_index_from_cdf( cdf, numeric::random::rg() );
	--bin_from_cdf;

	Size phi_ind = bin_from_cdf / n_psi_;
	Size psi_ind = bin_from_cdf - phi_ind * n_phi_;

	// following lines set phi and set to values drawn proportionately from Rama space
	// AS Nov 2013 - note that phi/psi bins are lower left corners, so we only want to add "noise" from a uniform distribution
	phi = binw_ * ( phi_ind + numeric::random::uniform() );
	psi = binw_ * ( psi_ind + numeric::random::uniform() );
}

} // scoring
} // core
