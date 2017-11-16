// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/RamaPrePro.cc
/// @brief
/// @author Frank DiMaio
/// @author Vikram K. Mulligan (vmullig@uw.edu) Feb 2016 -- made this compatible with canonical D-amino acids; returns 0 for noncanonicals.  Also, refactored greatly in Sept 2016 to support noncanonical rama tables with arbitrary numbers of degrees of freedom.

// Unit Headers
#include <core/scoring/RamaPrePro.hh>

// Package Headers
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/mainchain_potential/MainchainScoreTable.hh>
#include <basic/basic.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>
#include <numeric/MathMatrix.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <boost/algorithm/string.hpp>


namespace core {
namespace scoring {

static basic::Tracer TR("core.scoring.RamaPrePro");

RamaPrePro::RamaPrePro() :
	canonical_score_tables_(),
	canonical_prepro_score_tables_()
{
	using namespace basic::options;
	read_canonical_rpp_tables( );
}

/// @brief Evaluate the rama score for this residue (res1) given the identity of the next (res2).
/// @details This version only works for noncanonical or canonical residues with any number of mainchain
/// torsions.  If the next residue's identity is pro or d-pro, a different score table is used.  Note:
/// if return_derivs is true, the gradient vector is populated only.  If it is false, then only the
/// score_rama value is populated.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
RamaPrePro::eval_rpp_rama_score(
	core::conformation::Conformation const & conf,
	core::chemical::ResidueTypeCOP res1,
	core::chemical::ResidueTypeCOP res2,
	utility::vector1 < core::Real > mainchain_torsions, //Deliberately copied, not passed by reference
	Real & score_rama,
	utility::vector1 < core::Real > &gradient,
	bool const return_derivs
) const {
	core::chemical::AA const res_aa1( res1->backbone_aa() );
	if ( (core::chemical::is_canonical_L_aa(res_aa1) || core::chemical::is_canonical_D_aa(res_aa1) || res_aa1 == core::chemical::aa_gly ) &&
			!res1->defines_custom_rama_prepro_map( is_N_substituted( res2 ) )
			) {
		debug_assert( mainchain_torsions.size() == 2 );
		core::Real denergy_dphi, denergy_dpsi;
		eval_rpp_rama_score( res_aa1, res2, mainchain_torsions[1], mainchain_torsions[2], score_rama, denergy_dphi, denergy_dpsi, return_derivs ); //Call the version that uses fast enum-based lookups of scoring tables.
		if ( return_derivs ) {
			gradient.resize(2);
			gradient[1] = denergy_dphi; gradient[2] = denergy_dpsi;
		}
		return;
	}

	//Otherwise, we need to do slower scoring table lookups:
	bool const is_d( res1->is_d_aa() ); //Score tables are loaded for L-amino acids; D-versions are mirrored.

	//Flip torsions for D-amino acids.
	if ( is_d ) {
		for ( core::Size i=1, imax=mainchain_torsions.size(); i<=imax; ++i ) mainchain_torsions[i] *= -1.0;
	}

	core::chemical::ResidueTypeCOP ltype( is_d ? conf.residue_type_set_for_conf( res1->mode() )->get_mirrored_type( res1 )  : res1 );

	ScoringManager* manager( ScoringManager::get_instance() );

	core::chemical::mainchain_potential::MainchainScoreTableCOP cur_table( manager->get_rama_prepro_mainchain_torsion_potential( ltype, true, is_N_substituted( res2 ) ) );
	if ( !cur_table ) { //No scoring table defined for this residue type.
		if ( return_derivs ) {
			gradient.resize( res1->mainchain_atoms().size() - 1 );
			gradient[1] = 0.0;
			gradient[2] = 0.0;
		} else {
			score_rama = 0.0;
		}
		return;
	}

	// If we reach this point, a scoring table is defined.
	debug_assert( mainchain_torsions.size() == res1->mainchain_atoms().size() - 1 );
	if ( return_derivs ) {
		gradient.resize( res1->mainchain_atoms().size() - 1 );
		cur_table->gradient( mainchain_torsions, gradient );
		// Flip gradient for D-amino acids.
		if ( is_d ) {
			for ( core::Size i=1, imax=gradient.size(); i<=imax; ++i ) gradient[i] *= -1.0;
		}
	} else {
		score_rama = cur_table->energy( mainchain_torsions );
	}
}

/// @brief Evaluate the rama score for this residue (res_aa1) given the identity of the next (res2).
/// @details This version only works for canonical L-amino acids, canonical D-amino acids, or glycine.  If the next
/// residue's identity is pro or d-pro, a different score table is used.  Note:
/// if return_derivs is true, the gradient vector is populated only.  If it is false, then only the
/// score_rama value is populated.
/// @author Rewritten by Vikram K. Mulligan (vmullig@uw.edu).
void
RamaPrePro::eval_rpp_rama_score(
	AA const res_aa1,
	core::chemical::ResidueTypeCOP res2,
	Real const phi,
	Real const psi,
	Real & score_rama,
	Real & denergy_dphi,
	Real & denergy_dpsi,
	bool const return_derivs
) const {

	denergy_dphi = 0; denergy_dpsi = 0; // Needed to fix "may be used uninitialized" compiler errors for calling code.
	bool const is_d( core::chemical::is_canonical_D_aa(res_aa1) );
	core::Real const d_multiplier( is_d ? -1.0 : 1.0 );

	//If this is neither a canonical D-amino acid, nor a canonical L-amino acid, nor glycine return 0:
	if ( !core::chemical::is_canonical_L_aa( res_aa1 ) && !is_d && res_aa1 != core::chemical::aa_gly ) {
		if ( return_derivs ) {
			denergy_dphi = 0.0;
			denergy_dpsi = 0.0;
		} else {
			score_rama = 0.0;
		}
		return;
	}

	//Get the L-equivalent if this is a canonical D-residue:
	core::chemical::AA const res_aa1_copy( is_d ? core::chemical::get_L_equivalent(res_aa1) : res_aa1 );

	if ( res_aa1_copy > core::chemical::num_canonical_aas ) { //Noncanonical case: return 0.
		if ( return_derivs ) {
			denergy_dphi = 0.0;
			denergy_dpsi = 0.0;
		} else {
			score_rama = 0.0;
		}
	} else { //Canonical case: return something
		utility::vector1< core::Real > phipsi(2);
		phipsi[1] = d_multiplier * phi;
		phipsi[2] = d_multiplier * psi;
		utility::vector1 < core::Real > derivs(2);
		core::chemical::mainchain_potential::MainchainScoreTableCOP cur_table;
		if ( is_N_substituted( res2 ) ) { //VKM -- crude approximation: this residue is considered "pre-pro" if it precedes an L- or D-proline.  (The N and CD are achiral).
			cur_table =  canonical_prepro_score_tables_.at(res_aa1_copy);
		} else {
			cur_table =  canonical_score_tables_.at(res_aa1_copy);
		}
		debug_assert(cur_table);
		if ( return_derivs ) {
			cur_table->gradient( phipsi, derivs );
			denergy_dphi = d_multiplier * derivs[1];
			denergy_dpsi = d_multiplier * derivs[2];
		} else {
			score_rama = cur_table->energy(phipsi);
		}
	}
}

/// @brief Given the current residue (res1) and the next one (res2), randomly draw mainchain torsion values for the current
/// residue, biased by the Ramachandran probabilities for its type.
/// @details This version is general, usable for canonicals or noncanonicals.
/// @param[in] res1 The current residue, for which we're drawing mainchain torsions.
/// @param[in] res2 The next residue, used to determine whether to use pre-proline tables or not.
/// @param[out] torsions A vector of mainchain torsions for the current residue.
/// @note If the mainchain potential is a function of fewer than all of the mainchain torsions, the vector returned will
/// not have random values for the torsions that the mainchain potential is not a function of.  Rather, these torsions will
/// be set to 0.  Use the RamaPrePro::get_mainchain_torsions_covered() function to get a vector of torsion indices that were
/// randomized based on the mainchain potential CDF.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
RamaPrePro::random_mainchain_torsions(
	core::conformation::Conformation const & conf,
	core::chemical::ResidueTypeCOP res1,
	core::chemical::ResidueTypeCOP res2,
	utility::vector1 < core::Real > &torsions
) const {
	core::chemical::AA const res_aa1( res1->backbone_aa() );

	if ( (core::chemical::is_canonical_L_aa(res_aa1) || core::chemical::is_canonical_D_aa(res_aa1) || res_aa1 == core::chemical::aa_gly ) &&
			!res1->defines_custom_rama_prepro_map( is_N_substituted( res2 ) )
			) { //This is canonical.
		random_mainchain_torsions( res_aa1, res2, torsions );
		return;
	}

	//Otherwise, this is a noncanonical, and the logic becomes more complicated.
	//First, we must determine whether this is a D-amino acid.  (Score tables are for L-versions.)
	bool const is_d( res1->is_d_aa() );

	//Get the L-equivalent type:
	core::chemical::ResidueTypeCOP ltype( is_d ? conf.residue_type_set_for_conf( res1->mode() )->get_mirrored_type( res1 )  : res1 );

	//Get an instance of the ScoringManager, and the proper MainchainScoreTable:
	ScoringManager* manager( ScoringManager::get_instance() );
	core::chemical::mainchain_potential::MainchainScoreTableCOP cur_table(
		manager->get_rama_prepro_mainchain_torsion_potential( ltype, true, is_N_substituted( res2 ) )
	);
	runtime_assert_string_msg( cur_table, "Error in core::scoring::RamaPrePro::random_mainchain_torsions(): No mainchain score table for residue type " + res1->name() + " exists." );

	//Draw torsions:
	cur_table->draw_random_mainchain_torsion_values( torsions );

	//Correct chirality, if necessary:
	if ( is_d ) {
		for ( core::Size i=1, imax=torsions.size(); i<=imax; ++i ) {
			torsions[i] *= -1.0;
		}
	}
}

/// @brief Get a const reference to the vector of mainchain torsions indices that the mainchain potential for a given noncanonical ResidueType covers.
/// @details For example, for an oligourea, this would return {1, 2, 3}, since the Rama maps for oligoureas cover
/// phi, theta, and psi (mainchain torsions 1, 2, and 3, respectively), but not mu or omega (mainchain torsions 4 and 5).
/// @param[in] conf The current conformation, for reference.
/// @param[in] res1 The current residue, for which we're drawing mainchain torsions.
/// @param[in] res2 The next residue, used to determine whether to use pre-proline tables or not.
utility::vector1< core::Size > const &
RamaPrePro::get_mainchain_torsions_covered(
	core::conformation::Conformation const & conf,
	core::chemical::ResidueTypeCOP res1,
	core::chemical::ResidueTypeCOP res2
) const {
	//Otherwise, this is a noncanonical, and the logic becomes more complicated.
	//First, we must determine whether this is a D-amino acid.  (Score tables are for L-versions.)
	bool const is_d( res1->is_d_aa() );

	//Get the L-equivalent type:
	core::chemical::ResidueTypeCOP ltype( is_d ? conf.residue_type_set_for_conf( res1->mode() )->get_mirrored_type( res1 )  : res1 );

	//Get an instance of the ScoringManager, and the proper MainchainScoreTable:
	ScoringManager* manager( ScoringManager::get_instance() );
	core::chemical::mainchain_potential::MainchainScoreTableCOP cur_table(
		manager->get_rama_prepro_mainchain_torsion_potential( ltype, true, is_N_substituted( res2 ) )
	);
	runtime_assert_string_msg( cur_table, "Error in core::scoring::RamaPrePro::get_mainchain_torsions_covered(): No mainchain score table for residue type " + res1->name() + " exists." );

	return cur_table->get_mainchain_torsions_covered();
}

/// @brief Given the current residue (res1) and the next one (res2), randomly draw mainchain torsion values for the current
/// residue, biased by the Ramachandran probabilities for its type.
/// @details This version is for canonical residues only.
/// @param[in] res_aa1 The AA enum for a canonical residue type.
/// @param[in] res2 The next residue, used to determine whether to use pre-proline tables or not.
/// @param[out] torsions A vector of mainchain torsions for the current residue.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
RamaPrePro::random_mainchain_torsions(
	core::chemical::AA const res_aa1,
	core::chemical::ResidueTypeCOP res2,
	utility::vector1 < core::Real > &torsions
) const {

	//Determine whether this is a D-amino acid, and get the L-equivalent if it is:
	bool const is_d( core::chemical::is_canonical_D_aa(res_aa1) );
	core::chemical::AA const res_aa1_copy( is_d ? core::chemical::get_L_equivalent(res_aa1) : res_aa1 );
	runtime_assert_string_msg( res_aa1_copy <= core::chemical::num_canonical_aas,
		"Error in core::scoring::RamaPrePro::random_mainchain_torsions(): The canonical version of this function was called on a non-canonical amino acid." );

	//Get the appropriate MainchainScoreTable:
	core::chemical::mainchain_potential::MainchainScoreTableCOP cur_table;
	if ( is_N_substituted( res2 ) ) { //VKM -- crude approximation: this residue is considered "pre-pro" if it precedes an L- or D-proline.  (The N and CD are achiral).
		cur_table = canonical_prepro_score_tables_.at(res_aa1_copy);
	} else {
		cur_table = canonical_score_tables_.at(res_aa1_copy);
	}
	debug_assert(cur_table);

	//Call the random torsion generator:
	cur_table->draw_random_mainchain_torsion_values( torsions );

	//Correct chirality, if necessary:
	debug_assert( torsions.size() == 2 );
	if ( is_d ) {
		torsions[1] *= -1.0;
		torsions[2] *= -1.0;
	}
}

/// @brief Returns true if this aa is aa_pro or aa_dpr, false otherwise.
/// @details Also returns true for an N-methyl amino acid, peptoid, or
/// proline-like oligourea.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
bool
RamaPrePro::is_N_substituted(
	core::chemical::ResidueTypeCOP restype
) const {
	core::chemical::AA const aa( restype->aa() );
	return ( aa == core::chemical::aa_pro || aa == core::chemical::aa_dpr || aa == core::chemical::ou3_pro || restype->is_n_methylated() || restype->is_peptoid() );
}

/// @brief Ensure that the RamaPrePro scoring tables for the 20 canonical amino acids are set up, and that we are storing
/// pointers to them in a map of AA enum -> MainchainScoreTableCOP.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
RamaPrePro::read_canonical_rpp_tables( ) {

	canonical_score_tables_.clear();
	canonical_prepro_score_tables_.clear();

	ScoringManager* manager( ScoringManager::get_instance() ); //Raw pointer to the ScoringManager singleton.  Note that this is a rare instance in which there is no reason to use an owning pointer.
	core::chemical::ResidueTypeSetCOP rts( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

	//Read the score tables for this residue.
	for ( core::Size i( static_cast<core::Size>( core::chemical::first_l_aa )); i<=static_cast<core::Size>( core::chemical::num_canonical_aas ); ++i ) { //Loop through the canonical amino acids.
		core::chemical::ResidueTypeCOP restype( rts->get_representative_type_aa(static_cast< core::chemical::AA >(i) ) );
		runtime_assert( restype ); //There shouldn't be a problem getting the residue type.

		core::chemical::mainchain_potential::MainchainScoreTableCOP curtable( manager->get_rama_prepro_mainchain_torsion_potential( restype,  true /*RamaPrePro always uses polycubic interpolation*/, false /*Not a prepro table*/) );
		runtime_assert_string_msg( curtable, "Error in core::scoring::RamaPrePro::read_canonical_rpp_tables(): Could not read table for " + restype->name() + "." );
		canonical_score_tables_[ static_cast< core::chemical::AA >(i) ] = curtable;

		core::chemical::mainchain_potential::MainchainScoreTableCOP curtable_pp( manager->get_rama_prepro_mainchain_torsion_potential( restype,  true /*RamaPrePro always uses polycubic interpolation*/, true /*A prepro table*/ ) );
		runtime_assert_string_msg( curtable_pp, "Error in core::scoring::RamaPrePro::read_canonical_rpp_tables(): Could not read pre-proline table for " + restype->name() + "." );
		canonical_prepro_score_tables_[ static_cast< core::chemical::AA >(i) ] = curtable_pp;
	}
}

} //namespace scoring
} //namespace core
