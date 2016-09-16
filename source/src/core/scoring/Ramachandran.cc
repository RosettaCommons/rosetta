// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/Ramachandran.cc
/// @brief  Ramachandran potential class implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// @author Modified by Vikram K. Mulligan (vmullig@uw.edu) for D-amino acids, noncanonical alpha-amino acids, etc.
/// @details  This class should really be rewritten from the ground up, to permit greater flexibility and to remove protein-only assumptions that
/// limit the utility of the rama score term.  This is very old, kind of ugly, not very versatile code.

// Unit Headers
#include <core/scoring/Ramachandran.hh>

// Package Headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/ppo_torsion_bin.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/interpolation/periodic_range/half/interpolation.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>
#include <numeric/xyz.functions.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray4D.hh>
#include <ObjexxFCL/string.functions.hh>

// AS -- to get access to get_torsion_bin()
#include <core/conformation/util.hh>

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <utility/vector1.hh>
#include <sstream>

//MaximCode
#include <basic/Tracer.hh>
#include <boost/algorithm/string.hpp>

using namespace ObjexxFCL;

//MaximCode
static basic::Tracer TR("core.scoring.ramachandran");

namespace core {
namespace scoring {


// @brief Auto-generated virtual destructor
Ramachandran::~Ramachandran() = default;

typedef Ramachandran R;

Real const R::binw_( 10.0 );
Real const R::rama_sampling_thold_( 0.00075 ); // only sample torsions with Rama prob above this value
Real const R::rama_sampling_factor_( 10.0 ); // factor for increased precision of Rama sampling table


Ramachandran::Ramachandran() :
	ram_probabil_( n_phi_, n_psi_, 3, n_aa_ ),
	extra_ram_probabil_(),
	ram_counts_( n_phi_, n_psi_, 3, n_aa_ ),
	extra_ram_counts_(),
	ram_energ_( n_phi_, n_psi_, 3, n_aa_),
	extra_ram_energ_(),
	ram_entropy_( 3, n_aa_ ),
	cdf_( n_aa_ ),
	extra_cdf_(),
	cdf_by_torsion_bin_( n_aa_, conformation::n_ppo_torsion_bins ),
	n_valid_pp_bins_by_ppo_torbin_( n_aa_, conformation::n_ppo_torsion_bins ),
	phi_psi_bins_above_thold_( n_aa_ ),
	use_rama_power_(false),
	rama_power_(1.0)
{
	if ( basic::options::option[ basic::options::OptionKeys::score::rama_power ].user() ) {
		set_use_rama_power(true);
		set_rama_power( basic::options::option[ basic::options::OptionKeys::score::rama_power ]() );
	}
	if ( !basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable]
			|| !basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_rama_enable] ) {
		read_rama(
			basic::options::option[basic::options::OptionKeys::corrections::score::rama_map]().name(),
			basic::options::option[basic::options::OptionKeys::corrections::score::use_bicubic_interpolation]);
	} else {
		{
			std::string _smoothingRequsted = basic::options::option[ basic::options::OptionKeys::corrections::shapovalov_lib::shap_rama_smooth_level ];
			std::string _smoothingAsKeyword = "undefined";

			if ( _smoothingRequsted.compare("1")==0 ) {
				_smoothingAsKeyword = "lowest_smooth";
			} else if ( _smoothingRequsted.compare("2")==0 ) {
				_smoothingAsKeyword = "lower_smooth";
			} else if ( _smoothingRequsted.compare("3")==0 ) {
				_smoothingAsKeyword = "higher_smooth";
			} else if ( _smoothingRequsted.compare("4")==0 ) {
				_smoothingAsKeyword = "highest_smooth";
			} else {
				_smoothingAsKeyword = "unknown";
			}
			TR << "shapovalov_lib::shap_rama_smooth_level of " << _smoothingRequsted << "( aka " << _smoothingAsKeyword << " )" << " got activated." << std::endl;
		}

		read_rama(
			basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_rama_map]().name(),
			basic::options::option[basic::options::OptionKeys::corrections::score::use_bicubic_interpolation]);
	}
}

Ramachandran::Ramachandran(
	std::string const & rama_map_filename,
	bool use_bicubic_interpolation
) :
	ram_probabil_( n_phi_, n_psi_, 3, n_aa_ ),
	extra_ram_probabil_(),
	ram_counts_( n_phi_, n_psi_, 3, n_aa_ ),
	extra_ram_counts_(),
	ram_energ_( n_phi_, n_psi_, 3, n_aa_),
	extra_ram_energ_(),
	ram_entropy_( 3, n_aa_ ),
	cdf_( n_aa_ ),
	extra_cdf_(),
	cdf_by_torsion_bin_( n_aa_, conformation::n_ppo_torsion_bins ),
	n_valid_pp_bins_by_ppo_torbin_( n_aa_, conformation::n_ppo_torsion_bins ),
	phi_psi_bins_above_thold_( n_aa_ ),
	use_rama_power_(false),
	rama_power_(1.0)
{
	if ( basic::options::option[ basic::options::OptionKeys::score::rama_power ].user() ) {
		set_use_rama_power(true);
		set_rama_power( basic::options::option[ basic::options::OptionKeys::score::rama_power ]() );
	}
	read_rama(rama_map_filename, use_bicubic_interpolation);
}

/// @brief Given a custom Rama table type, get the filename of the database file containing
/// the data.
/// @details Accesses the options system for this information.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
std::string
Ramachandran::get_custom_rama_table_filename (
	Rama_Table_Type const type
) const  {
	using namespace basic::options;
	using namespace basic::options::OptionKeys::corrections::score;

	switch(type) {
	case(flat_l_aa_ramatable) :
		return option[rama_map_average_L_flat]();
	case(flat_d_aa_ramatable) :
		return option[rama_map_average_L_flat]();
	case(flat_symm_dl_aa_ramatable) :
		return option[rama_map_sym_average_L_flat]();
	case(flat_symm_gly_ramatable) :
		return option[rama_map_sym_gly_flat]();
	case(flat_symm_pro_ramatable) :
		return option[rama_map_sym_pro_flat]();
	case(flat_l_aa_ramatable_stringent) :
		return option[rama_map_average_L_flat_stringent]();
	case(flat_d_aa_ramatable_stringent) :
		return option[rama_map_average_L_flat_stringent]();
	case(flat_symm_dl_aa_ramatable_stringent) :
		return option[rama_map_sym_average_L_flat_stringent]();
	case(flat_symm_gly_ramatable_stringent) :
		return option[rama_map_sym_gly_flat_stringent]();
	case(flat_symm_pro_ramatable_stringent) :
		return option[rama_map_sym_pro_flat_stringent]();
	default :
		break;
	}

	utility_exit_with_message( "Error in core::scoring::Ramachandran::get_custom_rama_table_filename(): No filename associated with unknown rama table type." );
	return "";
}


/// @brief Given a Rama_Table_Type name, return the type.
/// @details Calls get_ramatable_name_by_type().
/// @author Vikram K. Mulligan (vmullig@uw.edu)
Rama_Table_Type
Ramachandran::get_ramatable_type_by_name (
	std::string const &name
) const {
	for ( core::Size i=1; i<=static_cast<core::Size>(end_of_ramatable_type_list); ++i ) {
		if ( get_ramatable_name_by_type( static_cast<Rama_Table_Type>(i) ) == name ) {
			return static_cast<Rama_Table_Type>(i);
		}
	}
	return unknown_ramatable_type;
}

/// @brief Given a Rama_Table_Type, return the name.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
std::string
Ramachandran::get_ramatable_name_by_type (
	Rama_Table_Type const type
) const {
	switch(type) {
	case(flat_l_aa_ramatable) :
		return "flat_l_aa_ramatable";
	case(flat_d_aa_ramatable) :
		return "flat_d_aa_ramatable";
	case(flat_symm_dl_aa_ramatable) :
		return "flat_symm_dl_aa_ramatable";
	case(flat_symm_gly_ramatable) :
		return "flat_symm_gly_ramatable";
	case(flat_symm_pro_ramatable) :
		return "flat_symm_pro_ramatable";
	case(flat_l_aa_ramatable_stringent) :
		return "flat_l_aa_ramatable_stringent";
	case(flat_d_aa_ramatable_stringent) :
		return "flat_d_aa_ramatable_stringent";
	case(flat_symm_dl_aa_ramatable_stringent) :
		return "flat_symm_dl_aa_ramatable_stringent";
	case(flat_symm_gly_ramatable_stringent) :
		return "flat_symm_gly_ramatable_stringent";
	case(flat_symm_pro_ramatable_stringent) :
		return "flat_symm_pro_ramatable_stringent";
	default :
		break;
	}
	return "unknown_ramatable_type";
}

///////////////////////////////////////////////////////////////////////////////

/// @brief Returns true if passed a core::chemical::AA corresponding to a
/// D-amino acid, and false otherwise.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
bool
Ramachandran::is_canonical_d_aminoacid(
	AA const res_aa
) const {
	return core::chemical::is_canonical_D_aa(res_aa);
}

///////////////////////////////////////////////////////////////////////////////

/// @brief When passed a d-amino acid, returns the l-equivalent.  Returns
/// aa_unk otherwise.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
core::chemical::AA
Ramachandran::get_l_equivalent(
	AA const d_aa
) const {
	return core::chemical::get_L_equivalent(d_aa);
}


///////////////////////////////////////////////////////////////////////////////

/// @brief evaluate rama score for each (protein) residue and store that score
/// in the pose.energies() object
void
Ramachandran::eval_rama_score_all(
	pose::Pose & pose,
	ScoreFunction const & scorefxn
) const
{
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

	int const total_residue = pose.size();

	// exclude chain breaks

	Energies & pose_energies( pose.energies() );

	for ( int ii = 1; ii <= total_residue; ++ii ) {
		if ( !pose.residue(ii).is_protein() || pose.residue(ii).is_terminus() || pose.residue(ii).is_virtual_residue() ) continue;

		Real rama_score,dphi,dpsi;
		if ( is_normally_connected(pose.residue(ii)) ) {
			eval_rama_score_residue(pose.residue(ii),rama_score,dphi,dpsi);
			//printf("Residue %i is normal.\n", ii); fflush(stdout); //DELETE ME -- FOR TESTING ONLY
			//std::cout << "Rama: residue " << ii << " = " << rama_score << std::endl;
			pose_energies.onebody_energies( ii )[rama] = rama_score;
		} else {
			//printf("Residue %i: THIS SHOULD HAPPEN ONLY IF THIS RESIDUE HAS WEIRD CONNECTIONS.", ii); fflush(stdout); //DELETE ME -- FOR TESTING ONLY
			eval_rama_score_residue_nonstandard_connection(pose, pose.residue(ii),rama_score,dphi,dpsi);
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
void
Ramachandran::write_rama_score_all( Pose const & /*pose*/ ) const
{}


///////////////////////////////////////////////////////////////////////////////
/// Sample phi/psi torsions with probabilities proportionate to their
/// Ramachandran probabilities
/// Note -- this function had previously required that the option
/// loops::nonpivot_torsion_sampling be active.  This function now
/// performs a just-in-time check to initialize these tables the first
/// time they are requested -- To properly multi-thread this code, the
/// function should nab a mutex so that no two threads try to execute
/// the code at once.
/// Note2 -- this function has been hackily patched to allow it to be used for
/// D-amino acids. (VKM -- 21 Aug 2014).
void
Ramachandran::random_phipsi_from_rama(
	AA const res_aa,
	Real & phi,
	Real & psi
) const {
	AA res_aa2 = res_aa;
	core::Real phipsi_multiplier=1.0;
	if ( is_canonical_d_aminoacid(res_aa) ) {
		res_aa2=get_l_equivalent(res_aa);
		phipsi_multiplier=-1.0;
	}
	draw_random_phi_psi_from_cdf( cdf_[ res_aa2 ], phi, psi );
	phi *= phipsi_multiplier;
	psi *= phipsi_multiplier;
	return;
}

void
Ramachandran::uniform_phipsi_from_allowed_rama(
	AA const res_aa,
	Real & phi,
	Real & psi
) const {
	using numeric::random::uniform;
	using numeric::random::random_range;

	Size n_torsions = phi_psi_bins_above_thold_[res_aa].size();
	Size random_index = random_range(1, n_torsions);
	Size ppbin_index = phi_psi_bins_above_thold_[ res_aa ][ random_index ];
	Size phi_index = ppbin_index / n_psi_;
	Size psi_index = ppbin_index - phi_index * n_psi_;

	phi = binw_ * ( phi_index + uniform());
	psi = binw_ * ( psi_index + uniform());
}

/// @brief return true if the bin that the given phi/psi pair falls in has a
/// probability that exceeds the rama_sampling_thold_ constant.
bool
Ramachandran::phipsi_in_allowed_rama(
	AA const aa,
	Real phi,
	Real psi
) const {
	phi = numeric::nonnegative_principal_angle_degrees(phi);
	psi = numeric::nonnegative_principal_angle_degrees(psi);

	Size phi_bin = static_cast< Size > ( phi / binw_ );
	Size psi_bin = static_cast< Size > ( psi / binw_ );

	return ram_probabil_(phi_bin+1, psi_bin+1, 3, aa ) > rama_sampling_thold_;
}

bool
Ramachandran::phipsi_in_forbidden_rama(
	AA const aa,
	Real phi,
	Real psi
) const {
	return ! phipsi_in_allowed_rama(aa, phi, psi);
}

///////////////////////////////////////////////////////////////////////////////
/// Sample phi/psi torsions with probabilities proportionate to their
/// Ramachandran probabilities -- this version performs lookup restricted to specified torsion bins
/// based on random_phipsi_from_rama and has the same issue for parallel running

/// @author Amelie Stein (amelie.stein@ucsf.edu)
/// @date Fri May 11 15:52:01 PDT 2012
/// @details returns a random phi/psi combination within the given torsion bin -- WARNING: this will only work for the torsion bins that are currently implemented

void
Ramachandran::random_phipsi_from_rama_by_torsion_bin(
	AA const res_aa,
	Real & phi,
	Real & psi,
	conformation::ppo_torsion_bin torsion_bin
) const {
	draw_random_phi_psi_from_cdf( cdf_by_torsion_bin_( res_aa, torsion_bin ), phi, psi );
}


void
Ramachandran::get_entries_per_torsion_bin(
	AA const res_aa,
	std::map< conformation::ppo_torsion_bin, core::Size > & tb_frequencies
) const {
	tb_frequencies[ conformation::ppo_torbin_A ] = phi_psi_bins_above_thold_[ res_aa ][ conformation::ppo_torbin_A ];
	tb_frequencies[ conformation::ppo_torbin_B ] = phi_psi_bins_above_thold_[ res_aa ][ conformation::ppo_torbin_B ];
	tb_frequencies[ conformation::ppo_torbin_E ] = phi_psi_bins_above_thold_[ res_aa ][ conformation::ppo_torbin_E ];
	tb_frequencies[ conformation::ppo_torbin_G ] = phi_psi_bins_above_thold_[ res_aa ][ conformation::ppo_torbin_G ];
	tb_frequencies[ conformation::ppo_torbin_X ] = phi_psi_bins_above_thold_[ res_aa ][ conformation::ppo_torbin_X ];
}


///////////////////////////////////////////////////////////////////////////////
void
Ramachandran::eval_rama_score_residue_nonstandard_connection(
	core::pose::Pose const & mypose,
	conformation::Residue const & res,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const {

	if ( res.connection_incomplete(1) || res.connection_incomplete(2) ) { //If this is an open terminus, don't score this residue.
		rama=0.0;
		drama_dphi=0.0;
		drama_dpsi=0.0;
		return;
	}

	core::conformation::Residue const & lowerres = mypose.residue(res.residue_connection_partner(1));
	core::conformation::Residue const & upperres = mypose.residue(res.residue_connection_partner(2));

	// NOTE: The following assumes that we're scoring something with an alpha-amino acid backbone!  At the time
	// of this writing (7 Feb 2013), rama checks that the residue being scored is either a standard L- or D-amino acid.
	Real phi = numeric::nonnegative_principal_angle_degrees(
		numeric::dihedral_degrees(
		lowerres.xyz( lowerres.residue_connect_atom_index(res.residue_connection_conn_id(1)) ),
		res.xyz("N"), //Position of N
		res.xyz("CA"), //Position of CA
		res.xyz("C") //Position of C
		));

	Real psi=numeric::nonnegative_principal_angle_degrees(
		numeric::dihedral_degrees(
		res.xyz("N"), //Position of N
		res.xyz("CA"), //Position of CA
		res.xyz("C"), //Position of C
		upperres.xyz( upperres.residue_connect_atom_index(res.residue_connection_conn_id(2)) )
		));

	//printf("rsd %lu phi=%.3f psi=%.3f\n", res.seqpos(), phi, psi); fflush(stdout); //DELETE ME

	core::chemical::AA ref_aa = res.aa();
	if ( res.backbone_aa() != core::chemical::aa_unk ) ref_aa = res.backbone_aa(); //If this is a noncanonical that specifies a canonical to use as a Rama template, use the template.

	eval_rama_score_residue( ref_aa, phi, psi, rama, drama_dphi, drama_dpsi );
}

///////////////////////////////////////////////////////////////////////////////
Real
Ramachandran::eval_rama_score_residue(
	conformation::Residue const & rsd
) const {
	Real rama, drama_dphi, drama_dpsi;
	eval_rama_score_residue( rsd, rama, drama_dphi, drama_dpsi );
	return rama;
}

///////////////////////////////////////////////////////////////////////////////
void
Ramachandran::eval_rama_score_residue(
	conformation::Residue const & rsd,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const {
	using namespace numeric;

	debug_assert( rsd.is_protein() );

	if ( 0.0 == nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1) ) ||
			0.0 == nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2) ) ||
			rsd.is_terminus() ||
			rsd.is_virtual_residue() ) { // begin or end of chain -- don't calculate rama score
		rama = 0.0;
		drama_dphi = 0.0;
		drama_dpsi = 0.0;
		return;
	}

	Real phi=0.0;
	Real psi=0.0;


	if ( is_normally_connected(rsd) ) { //If this residue is conventionally connected
		phi = nonnegative_principal_angle_degrees( rsd.mainchain_torsion(1) );
		psi = nonnegative_principal_angle_degrees( rsd.mainchain_torsion(2) );
		core::chemical::AA ref_aa = rsd.aa();
		if ( rsd.backbone_aa() != core::chemical::aa_unk ) {
			ref_aa = rsd.backbone_aa(); //If this is a noncanonical that specifies a canonical to use as a Rama template, use the template.
		}

		eval_rama_score_residue( ref_aa, phi, psi, rama, drama_dphi, drama_dpsi );
	} else { //If this residue is unconventionally connected (should be handled elsewhere)
		//printf("Residue %lu: THIS SHOULD NEVER OCCUR!\n", rsd.seqpos()); fflush(stdout); //DELETE ME -- FOR TESTING ONLY
		rama = 0.0;
		drama_dphi = 0.0;
		drama_dpsi = 0.0;
		return;
	}
}


///////////////////////////////////////////////////////////////////////////////
///
Real
Ramachandran::eval_rama_score_residue(
	AA const res_aa,
	Real const phi,
	Real const psi
) const {
	Real rama, drama_dphi, drama_dpsi;
	eval_rama_score_residue( res_aa, phi, psi, rama, drama_dphi, drama_dpsi );
	return rama;
}

void
Ramachandran::eval_rama_score_residue(
	AA const res_aa,
	Real const phi,
	Real const psi,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const {
	using namespace basic::options;
	eval_rama_score_residue(
		option[ OptionKeys::corrections::score::use_bicubic_interpolation ],
		option[ OptionKeys::corrections::score::rama_not_squared ],
		res_aa, phi, psi, rama, drama_dphi, drama_dpsi);
}

///////////////////////////////////////////////////////////////////////////////
///
void
Ramachandran::eval_rama_score_residue(
	bool use_bicubic_interpolation,
	bool rama_not_squared,
	AA const res_aa,
	Real const phi,
	Real const psi,
	Real & rama,
	Real & drama_dphi,
	Real & drama_dpsi
) const {
	using namespace numeric;

	core::chemical::AA res_aa2 = res_aa;
	core::Real phi2 = phi;
	core::Real psi2 = psi;
	core::Real d_multiplier=1.0; //A multiplier for derivatives: 1.0 for L-amino acids, -1.0 for D-amino acids.

	if ( is_canonical_d_aminoacid( res_aa ) ) { //If this is a D-amino acid, invert phi and psi and use the corresponding L-amino acid for the calculation
		res_aa2 = get_l_equivalent( res_aa );
		phi2 = -phi;
		psi2 = -psi;
		d_multiplier = -1.0;
	}

	if ( use_bicubic_interpolation ) {
		rama = rama_energy_splines_[ res_aa2 ].F(phi2,psi2);
		drama_dphi = d_multiplier*rama_energy_splines_[ res_aa2 ].dFdx(phi2,psi2);
		drama_dpsi = d_multiplier*rama_energy_splines_[ res_aa2 ].dFdy(phi2,psi2);

		if ( rama > 0.0 && use_rama_power() ) {
			//core::Real rama_power = basic::options::option[ basic::options::OptionKeys::score::rama_power ];
			// r = (r + 1) ** p - 1
			drama_dphi = rama_power() * pow(rama + 1.0, rama_power() - 1.0) * drama_dphi;
			drama_dpsi = rama_power() * pow(rama + 1.0, rama_power() - 1.0) * drama_dpsi;
			rama = pow(rama + 1.0, rama_power()) - 1.0;
		}

		return;
	} else {

		int ss_type = 3;

		FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
		FArray2A< Real > const rama_for_res( ram_probabil_(1, 1, ss_type, res_aa2), zero_index, zero_index );
		Real interp_p,dp_dphi,dp_dpsi;

		using namespace numeric::interpolation::periodic_range::half;
		interp_p = bilinearly_interpolated( phi2, psi2, binw_, n_phi_, rama_for_res, dp_dphi, dp_dpsi );

		if ( interp_p > 0.0 ) {
			rama = ram_entropy_(ss_type, res_aa2 ) - std::log( static_cast< double >( interp_p ) );
			double const interp_p_inv_neg = -1.0 / interp_p;
			drama_dphi = interp_p_inv_neg * d_multiplier * dp_dphi;
			drama_dpsi = interp_p_inv_neg * d_multiplier * dp_dpsi;
		} else {
			//if ( runlevel > silent ) { //apl fix this
			// std::cout << "rama prob = 0. in eval_rama_score_residue!" << std::endl;
			// std::cout << "phi" << SS( phi ) << " psi" << SS( psi ) <<
			//  " ss " << SS( ss ) << std::endl;
			//}
			drama_dphi = 0.0;
			drama_dpsi = 0.0;
			rama = 20.0;
		}

		if ( ! rama_not_squared ) {
			if ( rama > 1.0 ) {
				Real const rama_squared = rama * rama;
				if ( rama_squared > 20.0 ) {
					////  limit the score, but give the true derivative
					////   as guidance out of the flat section of map
					drama_dphi = 0.0;
					drama_dpsi = 0.0;
					rama = 20.0;
				} else {
					drama_dphi = 2 * rama * drama_dphi;
					drama_dpsi = 2 * rama * drama_dpsi;
					rama = rama_squared;
				}
			}
		}
	}
}

bool
Ramachandran::is_normally_connected (
	conformation::Residue const & res
) const {
	if ( res.is_upper_terminus() || res.is_lower_terminus() || res.connect_map_size() < 2 ) {
		return true; //Termini register as normally connected since they're not to be scored by rama.
	}

	bool firstconn = true, secondconn=true;

	if ( !res.connection_incomplete(1) ) firstconn= (res.residue_connection_partner(1) == res.seqpos()-1);
	if ( !res.connection_incomplete(2) ) secondconn= (res.residue_connection_partner(2) == res.seqpos()+1);

	return (firstconn && secondconn);
}

Size Ramachandran::n_phi_bins() const { return n_phi_; }

Size Ramachandran::n_psi_bins() const { return n_psi_; }

Real Ramachandran::rama_probability( core::chemical::AA aa, Real phi, Real psi ) const {
	phi = numeric::nonnegative_principal_angle_degrees( phi );
	psi = numeric::nonnegative_principal_angle_degrees( psi );
	Size phi_ind = static_cast< Size > ( floor( phi / binw_ ) ) + 1;
	Size psi_ind = static_cast< Size > ( floor( psi / binw_ ) ) + 1;
	return ram_probabil_( phi_ind, psi_ind, 3, aa );
}

Real Ramachandran::minimum_sampling_probability() const {
	return rama_sampling_thold_;
}


utility::vector1< Real > const &
Ramachandran::cdf_for_aa( core::chemical::AA aa ) const
{
	return cdf_[ aa ];
}

utility::vector1< Real > const &
Ramachandran::cdf_for_aa_for_torsion_bin(
	chemical::AA aa,
	conformation::ppo_torsion_bin torsion_bin
) const
{
	return cdf_by_torsion_bin_( aa, torsion_bin );
}

/// @brief Pick a random phi, psi value from a custom Rama table.
/// @details The custom Rama table is lazily loaded, so this function
/// is necessarily non-const.  By default, only the 20 canonical Rama
/// tables are loaded.
/// @param[in] type The type of custom rama table (an enum value).
/// @param[out] phi Randomly-drawn phi value, biased by the custom rama
/// table.
/// @param[out] psi Randomly-drawn psi value, biased by the custom rama
/// table.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
Ramachandran::draw_random_phi_psi_from_extra_cdf(
	Rama_Table_Type const type,
	Real & phi,
	Real & psi
) {
	bool is_d(false); //Are we dealing with a D-amino acid table?
	Rama_Table_Type type2(type);

	if ( type == flat_d_aa_ramatable ) {
		type2 = flat_l_aa_ramatable;
		is_d = true;
	} else if ( type == flat_d_aa_ramatable_stringent ) {
		type2 = flat_l_aa_ramatable_stringent;
		is_d = true;
	}

	if ( !has_custom_rama_cdf(type2) ) {
		generate_custom_rama_cdf(type2);
	}
	draw_random_phi_psi_from_cdf( custom_rama_cdf(type2), phi, psi );

	if ( is_d ) { //If this is a D-amino acid table, we've used the corresponding L-table, so we need to invert phi and psi.
		phi *= -1.0;
		psi *= -1.0;
	}
}

/// @brief If the -symmetric_gly_tables option is used, symmetrize the aa_gly table.
/// @details By default, the gly table is asymmetric because it is based on statistics from the PDB (which disproportionately put glycine
/// in the D-amino acid region of Ramachandran space).  However, the intrinsic propensities of glycine make it equally inclined to favour
/// right- or left-handed conformation.  (Glycine is achrial, and can't have a preference.)  Must be called AFTER gly table load, but prior
/// to bicubic interpolation setup.  Note that symmetrization is based on the probability table, and carries over to the energy table; the
/// counts table is left as-is (asymmetric).
/// @note If dont_use_shap is set to true, tables for all secondary structures are symmetrized.  If it's set to false, then only ss=3 tables
/// have been loaded, so they're the only ones symmetrized.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
Ramachandran::symmetrize_gly_table(
	bool const dont_use_shap
) {
	TR << "Symmetrizing glycine Ramachandran table." << std::endl;

	for ( core::Size ss=(dont_use_shap ? 1 : 3); ss<=3; ++ss ) { //Repeat for each secondary structure type, if we're not using the Shapovalov tables.  If we are, then only ss=3 has been loaded.

		//For debugging only:
		TR.Debug << "Old gly probability table for ss " << ss << ":" << std::endl;
		for ( core::Size iphi=1; iphi<=static_cast<core::Size>(n_phi_); ++iphi ) {
			for ( core::Size ipsi=1; ipsi<=static_cast<core::Size>(n_psi_); ++ipsi ) {
				TR.Debug << ram_probabil_(iphi,ipsi,ss,static_cast<int>(core::chemical::aa_gly)) << "\t";
			}
			TR.Debug << std::endl;
		}

		//Need to rebuild probability tables by averaging:
		core::Real probsum(0.0);
		for ( core::Size iphi=1; iphi<=static_cast<core::Size>(n_phi_); ++iphi ) {
			core::Size const opposite_phi( n_phi_ + 1 - iphi );
			for ( core::Size ipsi=1; ipsi<=static_cast<core::Size>(n_psi_); ++ipsi ) {
				core::Size const opposite_psi( n_psi_ + 1 - ipsi );
				core::Real const probval( (ram_probabil_(iphi,ipsi,ss,static_cast<int>(core::chemical::aa_gly)) + ram_probabil_(opposite_phi,opposite_psi,ss,static_cast<int>(core::chemical::aa_gly))) / 2.0 );
				ram_probabil_(iphi,ipsi,ss,static_cast<int>(core::chemical::aa_gly)) = probval;
				ram_probabil_(opposite_phi,opposite_psi,ss,static_cast<int>(core::chemical::aa_gly)) = probval;
				probsum += probval;
			}
		}

		//Need to loop through again to renormalize and calculate entropy:
		core::Real entropy(0.0);
		for ( core::Size iphi=1; iphi<=static_cast<core::Size>(n_phi_); ++iphi ) {
			for ( core::Size ipsi=1; ipsi<=static_cast<core::Size>(n_psi_); ++ipsi ) {
				ram_probabil_(iphi,ipsi,ss,static_cast<int>(core::chemical::aa_gly)) /= probsum;
				if ( ram_probabil_(iphi,ipsi,ss,static_cast<int>(core::chemical::aa_gly)) != 0 ) entropy += ram_probabil_(iphi,ipsi,ss,static_cast<int>(core::chemical::aa_gly)) * std::log( ram_probabil_(iphi,ipsi,ss,static_cast<int>(core::chemical::aa_gly)) );
			}
		}

		TR.Debug << "Old gly entropy for ss " << ss << " was " << ram_entropy_(ss,static_cast<int>(core::chemical::aa_gly)) << ".  New is " << entropy << "." << std::endl;

		ram_entropy_(ss,static_cast<int>(core::chemical::aa_gly)) = entropy;

		//For debugging only:
		TR.Debug << "New gly probability table for ss " << ss << ":" << std::endl;
		for ( core::Size iphi=1; iphi<=static_cast<core::Size>(n_phi_); ++iphi ) {
			for ( core::Size ipsi=1; ipsi<=static_cast<core::Size>(n_psi_); ++ipsi ) {
				TR.Debug << ram_probabil_(iphi,ipsi,ss,static_cast<int>(core::chemical::aa_gly)) << "\t";
			}
			TR.Debug << std::endl;
		}
		TR.Debug << "Old gly energy table for ss " << ss << ":" << std::endl;
		for ( core::Size iphi=1; iphi<=static_cast<core::Size>(n_phi_); ++iphi ) {
			for ( core::Size ipsi=1; ipsi<=static_cast<core::Size>(n_psi_); ++ipsi ) {
				TR.Debug << ram_energ_(iphi,ipsi,ss,static_cast<int>(core::chemical::aa_gly)) << "\t";
			}
			TR.Debug << std::endl;
		}

		//Finally, need to recalculate energy:
		for ( core::Size iphi=1; iphi<=static_cast<core::Size>(n_phi_); ++iphi ) {
			for ( core::Size ipsi=1; ipsi<=static_cast<core::Size>(n_psi_); ++ipsi ) {
				core::Real const pval( ram_probabil_(iphi,ipsi,ss,static_cast<int>(core::chemical::aa_gly)) );
				if ( pval != 0 ) {
					ram_energ_(iphi,ipsi,ss,static_cast<int>(core::chemical::aa_gly)) = -1.0*std::log( pval ) + entropy;
				} else {
					ram_energ_(iphi,ipsi,ss,static_cast<int>(core::chemical::aa_gly)) = 20.0; //Fixing this at 20 for the prohibited region for now.
				}
			}
		}

		TR.Debug << "New gly energy table for ss " << ss << ":" << std::endl;
		for ( core::Size iphi=1; iphi<=static_cast<core::Size>(n_phi_); ++iphi ) {
			for ( core::Size ipsi=1; ipsi<=static_cast<core::Size>(n_psi_); ++ipsi ) {
				TR.Debug << ram_energ_(iphi,ipsi,ss,static_cast<int>(core::chemical::aa_gly)) << "\t";
			}
			TR.Debug << std::endl;
		}

	} // Looping through secondary structures

	if ( TR.visible() ) { TR.flush(); }
	if ( TR.Debug.visible() ) { TR.Debug.flush(); }

	return;
}

void
Ramachandran::read_rama_map_file (
	utility::io::izstream * iunit
) {
	int aa_num,phi_bin,psi_bin,ss_type;
	Real check,min_prob,max_prob;
	double entropy;
	char line[60];
	int scan_count;
	float pval, eval; // vars for sscanf float I/O

	for ( int i = 1; i <= n_aa_ ; ++i ) {
		for ( int ii = 1; ii <= 3; ++ii ) {
			entropy = 0.0;
			check = 0.0;
			min_prob = 1e36;
			max_prob = -min_prob;
			for ( int j = 1; j <= 36; ++j ) {
				for ( int k = 1; k <= 36; ++k ) {
					iunit->getline( line, 60 );
					if ( iunit->eof() ) {
						return;
					} else if ( iunit->fail() ) { // Clear and continue: NO ERROR DETECTION
						iunit->clear();
					}
					std::sscanf( line, "%5d", &aa_num );
					std::sscanf( line+6, "%5d", &ss_type );
					std::sscanf( line+12, "%5d", &phi_bin );
					std::sscanf( line+18, "%5d", &psi_bin );
					std::sscanf( line+24, "%5d", &ram_counts_(j,k,ii,i) );
					std::sscanf( line+30, "%12f", &pval );
					ram_probabil_(j,k,ii,i) = pval;
					scan_count = std::sscanf( line+43, "%12f", &eval );
					ram_energ_(j,k,ii,i) = eval;

					if ( scan_count == EOF ) continue; // Read problem: NO ERROR DETECTION

					// This is the Slick & Slow (S&S) stream-based method that is too slow for large
					// files like this one, at least under the GCC 3.3.1 stream implementation.
					// It should be retried on future releases and target compilers because there is
					// no reason it cannot be competitive with good optimization and inlining.
					// If this is used the <cstdio> can be removed.
					//
					//     iunit >> bite( 5, aa_num ) >> skip( 1 ) >>
					//      bite( 5, ss_type ) >> skip( 1 ) >>
					//      bite( 5, phi_bin ) >> skip( 1 ) >>
					//      bite( 5, psi_bin ) >> skip( 1 ) >>
					//      bite( 5, ram_counts(j,k,ii,i) ) >> skip( 1 ) >>
					//      bite( 12, ram_probabil(j,k,ii,i) ) >> skip( 1 ) >>
					//      bite( 12, ram_energ(j,k,ii,i) ) >> skip;
					//     if ( iunit.eof() ) {
					//      goto L100;
					//     } else if ( iunit.fail() ) { // Clear and continue: NO ERROR DETECTION
					//      iunit.clear();
					//      iunit >> skip;
					//     }

					check += ram_probabil_(j,k,ii,i);
					entropy += ram_probabil_(j,k,ii,i) *
						std::log( static_cast< double >( ram_probabil_(j,k,ii,i) ) );
					min_prob = std::min(ram_probabil_(j,k,ii,i),min_prob);
					max_prob = std::max(ram_probabil_(j,k,ii,i),max_prob);
				}
			}
			ram_entropy_(ii,i) = entropy;
		}
		//cj  std::cout << SS( check ) << SS( std::log(min_prob) ) <<
		//cj   SS( std::log(max_prob) ) << SS( entropy ) << std::endl;
	}
}

void
Ramachandran::read_rama(
	std::string const & rama_map_filename,
	bool use_bicubic_interpolation
) {
	utility::io::izstream  iunit;

	// search in the local directory first
	iunit.open( rama_map_filename );

	if ( !iunit.good() ) {
		iunit.close();
		if ( !basic::database::open( iunit, rama_map_filename ) ) {
			std::stringstream err_msg;
			err_msg << "Unable to open Ramachandran map '" << rama_map_filename << "'.";
			utility_exit_with_message(err_msg.str());
		}
	}

	//cj      std::cout << "index" << "aa" << "ramachandran entropy" << std::endl;
	//KMa add_phospho_ser 2006-01

	bool const dont_use_shap(
		!basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable]
		|| !basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib::shap_rama_enable]
	);
	bool const symm_gly( basic::options::option[ basic::options::OptionKeys::score::symmetric_gly_tables ]() );

	if ( dont_use_shap ) {
		read_rama_map_file(&iunit);
	} else {
		read_rama_map_file_shapovalov(&iunit);
	}

	iunit.close();
	iunit.clear();

	//If the option is set to do this, symmetrize the glycine Ramachandran map.
	if ( symm_gly ) symmetrize_gly_table( dont_use_shap );

	if ( use_bicubic_interpolation ) {
		using namespace numeric;
		using namespace numeric::interpolation::spline;
		rama_energy_splines_.resize( chemical::num_canonical_aas );
		for ( Size ii = 1; ii <= chemical::num_canonical_aas; ++ii ) {
			BicubicSpline ramaEspline;
			MathMatrix< Real > energy_vals( 36, 36 );
			for ( Size jj = 0; jj < 36; ++jj ) {
				for ( Size kk = 0; kk < 36; ++kk ) {
					//MaximCode:
					if ( true ) {
						if ( ram_probabil_(jj+1,kk+1,3,ii ) != 0 ) {
							energy_vals( jj, kk ) = -std::log( ram_probabil_(jj+1,kk+1,3,ii )) + ram_entropy_(3,ii) ;
						} else {
							energy_vals( jj, kk ) = 20.0; //Arbitrarily setting the energy to +20 for cases of probability 0.
						}
					} else {
						energy_vals( jj, kk ) = ram_probabil_(jj+1,kk+1,3,ii );
					}
				}
			}
			BorderFlag periodic_boundary[2] = { e_Periodic, e_Periodic };
			Real start_vals[2];
			if (
					(
					basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable] &&
					basic::options::option[ basic::options::OptionKeys::corrections::shapovalov_lib::shap_rama_nogridshift ]
					)
					||
					(
					basic::options::option[basic::options::OptionKeys::corrections::shapovalov_lib_fixes_enable] &&
					basic::options::option[basic::options::OptionKeys::corrections::correct]
					)
					) {
				start_vals[0] = start_vals[1] = 0.0; // if this flag is on, then the Rama table is not shifted and aligns with the ten-degree boundaries.
			} else {
				start_vals[0] = start_vals[1] = 5.0; // otherwise, the grid is shifted by five degrees.
			}

			// Special case: don't shift the gly tables if the -symmetric_gly_tables flag is passed.
			// (This is not strictly correct if the Shapovalov tables are used, but good enough.)
			if ( ii == static_cast<core::Size>( core::chemical::aa_gly ) && symm_gly ) { start_vals[0] = start_vals[1] = 5.0; }

			Real deltas[2] = {10.0, 10.0}; // grid is 10 degrees wide
			bool lincont[2] = {false,false}; //meaningless argument for a bicubic spline with periodic boundary conditions
			std::pair< Real, Real > unused[2];
			unused[0] = std::make_pair( 0.0, 0.0 );
			unused[1] = std::make_pair( 0.0, 0.0 );
			ramaEspline.train( periodic_boundary, start_vals, deltas, energy_vals, lincont, unused );
			rama_energy_splines_[ ii ] = ramaEspline;
		}
	}

	initialize_rama_sampling_tables();
}

/// @brief Load a custom Ramachandran table, in addition to the 20x3 standard
/// ones that are always loaded.
/// @details Intended for sampling with alternative Ramachandran distributions.  Custom
/// tables are lazily loaded so as not to add to total Rosetta memory footprint.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
Ramachandran::load_custom_rama_table(
	Rama_Table_Type const type
) {
	runtime_assert_string_msg( !has_custom_rama_energy_table(type) && !has_custom_rama_probability_table(type) && !has_custom_rama_count_table(type),
		"Error in core::scoring::Ramachandran::load_custom_rama_table(): the tables for type " + get_ramatable_name_by_type(type) + " have already been loaded!" );

	std::string const filename( get_custom_rama_table_filename( type ) );

	utility::io::izstream iunit;

	// search in the local directory first
	iunit.open( filename );

	if ( !iunit.good() ) {
		iunit.close();
		if ( !basic::database::open( iunit, filename ) ) {
			std::stringstream err_msg("");
			err_msg << "Unable to open custom Ramachandran map \"" << filename << "\".";
			utility_exit_with_message(err_msg.str());
		}
	}

	if ( TR.visible() ) {
		TR << "Reading custom Ramachandran table from " << filename << "." << std::endl;
		TR.flush();
	}

	// Read/parse buffers:
	core::Size aa_num(0), phi_bin(0), psi_bin(0), ss_type (0);
	core::Real probsum(0.0);
	std::string line;

	//Storage for data that will be stored in extra_ram_probabil_, extra_ram_counts_, and extra_ram_energ_.
	utility::vector1 < utility::vector1 <core::Real> > prob;
	utility::vector1 < utility::vector1 <core::Size> > counts;
	utility::vector1 < utility::vector1 <core::Real> > energy;
	for ( core::Size j=1; j<=static_cast<core::Size>(n_phi_); ++j ) {
		utility::vector1 < core::Real > tempvect( n_psi_, 0.0 );
		utility::vector1 < core::Size > tempvect_size( n_psi_, 0.0 );
		prob.push_back( tempvect );
		counts.push_back( tempvect_size );
		energy.push_back( tempvect );
	}

	for ( core::Size j = 1; j <= static_cast<core::Size>(n_phi_); ++j ) {
		for ( core::Size k = 1; k <= static_cast<core::Size>(n_psi_); ++k ) {
			iunit.getline( line );
			if ( iunit.eof() ) {
				return;
			} else if ( iunit.fail() ) {
				std::stringstream err_msg("");
				err_msg << "Error reading file \"" << filename << "\".";
				utility_exit_with_message( err_msg.str() );
			}
			std::stringstream l(line);

			l >> aa_num;
			l >> ss_type;
			l >> phi_bin;
			l >> psi_bin;
			l >> counts[j][k];
			l >> prob[j][k];
			l >> energy[j][k];
			probsum += prob[j][k];

		}
	}
	iunit.close();

	for ( core::Size j=1; j<=static_cast<core::Size>(n_phi_); ++j ) {
		for ( core::Size k=1; k<=static_cast<core::Size>(n_psi_); ++k ) {
			prob[j][k] /= probsum; //Normalize the probabilities.
		}
	}

	// Store the loaded values:
	extra_ram_probabil_[type] = prob;
	extra_ram_counts_[type] = counts;
	extra_ram_energ_[type] = energy;
}

void
Ramachandran::initialize_rama_sampling_tables()
{
	init_rama_sampling_table( conformation::ppo_torbin_X );
	init_rama_sampling_table( conformation::ppo_torbin_A );
	init_rama_sampling_table( conformation::ppo_torbin_B );
	init_rama_sampling_table( conformation::ppo_torbin_E );
	init_rama_sampling_table( conformation::ppo_torbin_G );

	init_uniform_sampling_table();
}

void
Ramachandran::init_rama_sampling_table( conformation::ppo_torsion_bin torsion_bin )
{
	utility::vector1< Real > inner_cdf( n_phi_ * n_psi_, 0.0 );
	FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
	for ( int ii = 1; ii <= n_aa_; ++ii ) {
		std::fill( inner_cdf.begin(), inner_cdf.end(), Real( 0.0 ) );
		FArray2A< Real > const ii_rama_prob( ram_probabil_( 1, 1, 3, ii ), zero_index, zero_index );

		Size actual_allowed( 0 );
		Real allowed_probability_sum( 0.0 );
		for ( int jj = 0; jj < n_phi_; ++jj ) {
			for ( int kk = 0; kk < n_psi_; ++kk ) {

				inner_cdf[ jj*n_psi_ + kk +  1 ] = allowed_probability_sum;

				Real jjkk_prob = ii_rama_prob( jj, kk );
				if ( jjkk_prob < rama_sampling_thold_ ) continue;
				Real const cur_phi = binw_ * ( jj - ( jj > n_phi_ / 2 ? n_phi_ : 0 ));
				Real const cur_psi = binw_ * ( kk - ( kk > n_psi_ / 2 ? n_psi_ : 0 ));

				conformation::ppo_torsion_bin cur_tb = conformation::ppo_torbin_X;
				if ( torsion_bin != conformation::ppo_torbin_X ) {
					//  AS -- how can we get the factor properly / without hard-coding? - also: this takes very long...
					cur_tb = core::conformation::get_torsion_bin(cur_phi, cur_psi);
				}

				//std::cout << "    " << ii << " " << cur_phi << " " << cur_psi << " " << core::conformation::map_torsion_bin_to_char( cur_tb ) << " " << jjkk_prob << std::endl;

				if ( cur_tb == torsion_bin ) {
					++actual_allowed;
					// store the cumulative sum up to the current table entry
					// then increment the cumulative sum with the current table entry
					allowed_probability_sum += jjkk_prob;
				}
			}
		}

		//std::cout << "Ramachandran: " << ii << " " << core::conformation::map_torsion_bin_to_char( torsion_bin ) << " " << actual_allowed << std::endl;

		if ( actual_allowed == 0 ) {
			// store a vector of 0s instead of creating a vector of NaNs.
			cdf_by_torsion_bin_( ii, torsion_bin ) = inner_cdf;
			continue;
		}

		// now normalize the cdf vector by scaling by 1/allowed_probability_sum so that it represents a CDF that sums to 1.
		// and mark all zero-probability bins with the CDF value of their predecessor.
		Real const inv_allowed_probability_sum = 1 / allowed_probability_sum;
		for ( int jj = 1; jj <= n_phi_ * n_psi_; ++jj ) {
			inner_cdf[ jj ] *= inv_allowed_probability_sum;
		}

		// now update the cdf, cdf_by_torsion_bin, and n_valid_pp_bins_by_ppo_torbin tables
		if ( torsion_bin == conformation::ppo_torbin_X ) {
			cdf_[ ii ] = inner_cdf;
		}
		cdf_by_torsion_bin_( ii, torsion_bin ) = inner_cdf;
		n_valid_pp_bins_by_ppo_torbin_( ii, torsion_bin ) = actual_allowed;
	}
}

void
Ramachandran::init_uniform_sampling_table()
{
	FArray2A< Real >::IR const zero_index( 0, n_phi_ - 1);
	for ( int ii = 1; ii <= n_aa_; ++ii ) {
		utility::vector1< Size > minimum_prob_bin_inds;
		minimum_prob_bin_inds.reserve( n_phi_ * n_psi_ );
		FArray2A< Real > const ii_rama_prob( ram_probabil_( 1, 1, 3, ii ), zero_index, zero_index );
		for ( int jj = 0; jj < n_phi_; ++jj ) {
			for ( int kk = 0; kk < n_psi_; ++kk ) {
				Real jjkk_prob = ii_rama_prob( jj, kk );
				if ( jjkk_prob >= rama_sampling_thold_ ) {
					minimum_prob_bin_inds.push_back( jj*n_psi_ + kk );
				}
			}
		}
		phi_psi_bins_above_thold_[ ii ] = minimum_prob_bin_inds;
	}
}

/// @brief Pick a random phi and psi value given a cumulative distribution function.
/// @param[in] cdf The cumulative distribution function.
/// @param[out] phi The output phi value.
/// @param[out] psi The output psi value.
/// @details A random bin from the cumulative distribution function is chosen first.  Then uniform
/// randomness is added within the chosen bin to produce the final phi and psi values.
void
Ramachandran::draw_random_phi_psi_from_cdf(
	utility::vector1< Real > const & cdf,
	Real & phi,
	Real & psi
) const
{
	// the bin index can be unpacked to give the phi and psi indices
	Size bin_from_cdf = numeric::random::pick_random_index_from_cdf( cdf, numeric::random::rg() );
	--bin_from_cdf;

	Size phi_ind = bin_from_cdf / n_psi_;
	Size psi_ind = bin_from_cdf - phi_ind * n_phi_;

	//std::cout << " phi_ind: " << phi_ind << " psi_ind: " << psi_ind << std::endl;

	// following lines set phi and set to values drawn proportionately from Rama space
	// AS Nov 2013 - note that phi/psi bins are lower left corners, so we only want to add "noise" from a uniform distribution
	phi = binw_ * ( phi_ind + numeric::random::uniform() );
	psi = binw_ * ( psi_ind + numeric::random::uniform() );
}

/// @brief Generate a custom Rama cumulative distribution function from the corresponding energy table.
/// @details This is generated from the PROBABILITY table, not from the ENERGY table.  If the energy table has not been
/// loaded, this function loads the energy table first.
void
Ramachandran::generate_custom_rama_cdf(
	Rama_Table_Type const type
) {
	//Check that type refers to something sensible:
	runtime_assert_string_msg( type >= 1 && type < unknown_ramatable_type, "Error in core::scoring::Ramachandran::generate_custom_rama_cdf(): The type of Ramachandran table is not in the known types list." );

	//If we haven't loaded the energy table yet, load it now:
	if ( !has_custom_rama_probability_table( type ) ) {
		load_custom_rama_table( type );
	}

	//Compute the cumulative distribution function (inner_cdf):
	utility::vector1< Real > inner_cdf( n_phi_ * n_psi_, 0.0 );
	core::Real probsum(0.0);
	for ( core::Size jj = 0; jj < static_cast<core::Size>(n_phi_); ++jj ) {
		for ( core::Size kk = 0; kk < static_cast<core::Size>(n_psi_); ++kk ) {
			inner_cdf[ jj*n_psi_ + kk +  1 ] = probsum;
			probsum += custom_rama_probability_table(type, jj, kk);
		}
	}
	for ( core::Size i=1, imax=inner_cdf.size(); i<=imax; ++i ) {
		inner_cdf[i] = inner_cdf[i] / probsum;
	}

	//Store the computed CDF in the map:
	extra_cdf_[type] = inner_cdf;
}

//MaximCode
void
Ramachandran::read_rama_map_file_shapovalov (
	utility::io::izstream * iunit
) {
	//...
	//# Analysis      dataPoints-adaptive KDE
	//#
	//# res   phi     psi     Probabil        -log(Probabil)
	//#
	//ala     -180.0  -180.0  5.436719e-004   7.517165e+000
	//ala     -170.0  -180.0  9.476649e-004   6.961510e+000
	//ala     -160.0  -180.0  1.274955e-003   6.664844e+000
	//ala     -150.0  -180.0  1.334395e-003   6.619277e+000

	char line[256];
	Real check, min_prob, max_prob;

	int i; // aa index
	int ii = 3; // secondary structure index
	int j, k; //phi and psi indices
	char aa[4];
	std::string aaStr;
	double phi, psi, prob, minusLogProb;

	do {
		iunit->getline( line, 255 );

		if ( iunit->eof() ) {
			break;
		}

		if ( line[0]=='#' || line[0]=='\n' || line[0]=='\r' ) {
			continue;
		}

		std::sscanf( line, "%3s%lf%lf%lf%lf", aa, &phi, &psi, &prob, &minusLogProb );
		std::string prevAAStr = aaStr;
		aaStr = std::string(aa);
		boost::to_upper(aaStr);
		i = core::chemical::aa_from_name(aaStr);
		if ( aaStr != prevAAStr ) {
			if ( false ) {
				TR << "MaximCode test: " << "aa " << prevAAStr << "\tsumP " << check << "\tminP " << min_prob << "\tmaxP " << max_prob << std::endl;
			}

			check = 0.0;
			min_prob = 1e36;
			max_prob = -min_prob;
			ram_entropy_(ii,i) = 0;
		}

		//j, k (aka phi and psi indices) start with 1 and correspond to 0 and go all way to 36 corresponding to 350
		//since the storage type is FArray4D, Fortran-Compatible 4D Array where indices start from 1 by default
		if ( phi < 0 ) phi += 360;
		if ( psi < 0 ) psi += 360;

		//round produces unidentified symbol on Windows compilation, so replaced with ceil
		j = (int) ceil(phi / 10.0 - 0.5) + 1;
		k = (int) ceil(psi / 10.0 - 0.5) + 1;

		ram_probabil_(j,k,ii,i) = prob;
		ram_energ_(j,k,ii,i) = minusLogProb;
		//Maxim Shapovalov:
		//It is not entropy in accordance with textbooks.
		//It should be reversed with additional minus.
		//There is additional minus G = H - T*S in free energy
		//but it is skipped in the code by previous developers
		//So it means good.
		ram_entropy_(ii,i) += ram_probabil_(j,k,ii,i) * -ram_energ_(j,k,ii,i);

		check += ram_probabil_(j,k,ii,i);
		min_prob = std::min(ram_probabil_(j,k,ii,i),min_prob);
		max_prob = std::max(ram_probabil_(j,k,ii,i),max_prob);

	} while (true);
}


}
}
