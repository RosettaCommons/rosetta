// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/init/score_function_corrections.hh
/// @brief  Initialize score function corrections
/// @author Matthew O'meara


#ifndef INCLUDED_core_init_score_function_corrections_hh
#define INCLUDED_core_init_score_function_corrections_hh

#include <utility/options/OptionCollection.fwd.hh>
#include <string>

namespace core {
namespace init {

/// @brief Initialize the revert-to-Talaris2014
void init_revert_to_talaris( utility::options::OptionCollection & options );

/// @brief Reset a set of flags to their pre-talaris behavior
void revert_to_pre_talaris_2013_defaults( utility::options::OptionCollection & options );

/// @brief Initialize the "revert to before Talaris2013" mistake
void init_revert_to_pre_talaris_2013_mistake( utility::options::OptionCollection & options );

/// @brief Initialize the hbond Sp2 correction
/// Deprecated as Talaris2013 becomes default.
void init_hbond_sp2_correction( utility::options::OptionCollection & options );

/// @brief Initialize the FACTS correction
void init_facts_correction( utility::options::OptionCollection & options );

/// @brief Initialize the -correct correction
void init_correct_correction( utility::options::OptionCollection & options );

/// @brief Initialize the -shapovalov_lib_fixes_enable correction
void init_shapovalov_lib_fixes_enable_correction( utility::options::OptionCollection & options );

/// @brief Initialize the _most recent_ beta score function (-beta)
void init_beta_correction( utility::options::OptionCollection & options );

/// @brief Initialize the -beta_july15 score function
void init_beta_july15_correction( utility::options::OptionCollection & options );

/// @brief Initialize the -beta_nov15 score function; no longer valid
//void init_beta_nov15_correction( utility::options::OptionCollection & options );

/// @brief Initialize the -beta_nov16 score function
void init_beta_nov16_correction( utility::options::OptionCollection & options );

/// @brief Initialize the crystal refinement correction
void init_crystal_refinement_correction( utility::options::OptionCollection & options );

/// @brief Initialize nonideal bond geometry correction
void init_nonideal_correction( utility::options::OptionCollection & options );

/// @brief restore the the score function to the Score12prime version
void init_restore_score12prime( utility::options::OptionCollection & options );

/// @brief Alter the water/water and water/non-water solvation potential if the score:water_hybrid_sf flag is
/// on the command line
void init_spades_score_function_correction( utility::options::OptionCollection & options );

/// @brief Initialize the latest and greatest score function parameters
void init_score_function_corrections( utility::options::OptionCollection & options );

/// @brief Check if a score function is requested with incompatible option flags
/// Will return true if scorefunction is "sane" and false if not.
/// If throw_exception is true, will raise an exception instead of returning false.
bool check_score_function_sanity(
	utility::options::OptionCollection const & options,
	std::string const & scorefxn_key,
	bool throw_exception = false );

/// @brief  Apply some DNA-specific mods that are still in testing phase; only if -corrections::newdna present
void
init_dna_correction( utility::options::OptionCollection & options );


} // namespace
} // namespace

#endif // include_guard
