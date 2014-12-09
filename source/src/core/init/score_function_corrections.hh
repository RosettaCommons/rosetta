// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/init/score_function_corrections.hh
/// @brief  Initialize score function corrections
/// @author Matthew O'meara
///


#ifndef INCLUDED_core_init_score_function_corrections_hh
#define INCLUDED_core_init_score_function_corrections_hh

#include <string>

namespace core {
namespace init {

/// @brief Reset a set of flags to their pre-talaris behavior
void revert_to_pre_talaris_2013_defaults();

/// @brief Initialize the "revert to before Talaris2013" mistake
void init_revert_to_pre_talaris_2013_mistake();

/// @brief Initialize the hbond Sp2 correction
/// Deprecated as Talaris2013 becomes default.
void init_hbond_sp2_correction();

/// @brief Initialize the FACTS correction
void init_facts_correction();

/// @brief Initialize the -correct correction
void init_correct_correction();

/// @brief Initialize the -beta score function
void init_beta_correction();

/// @brief Initialize the crystal refinement correction
void init_crystal_refinement_correction();

/// @brief Initialize nonideal bond geometry correction
void init_nonideal_correction();

/// @brief restore the the score function to the Score12prime version
void init_restore_score12prime();

/// @brief Initialize the latest and greatest score function parameters
void init_score_function_corrections();

/// @brief Check if a score function is requested with incompatible option flags
void check_score_function_sanity(
	std::string const & scorefxn_key,
	bool warn_only=false);

/// @brief  Apply some DNA-specific mods that are still in testing phase; only if -corrections::newdna present
void
init_dna_correction();


} // namespace
} // namespace

#endif // include_guard
