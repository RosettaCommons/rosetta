// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/options/keys/OptionKeys.hh
/// @brief  utility::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note   Sample/starter OptionKey collection: Copy this file
///         to your project and adapt to the project namespace
///         and add the project OptionKeys


#ifndef INCLUDED_utility_options_keys_OptionKeys_hh
#define INCLUDED_utility_options_keys_OptionKeys_hh


// Package headers
#include <utility/options/keys/BooleanOptionKey.hh>


namespace utility {
namespace options {
namespace OptionKeys {


/// @brief Help option keys
extern BooleanOptionKey const help; // Show program help and exit


/// @brief Option display option keys
namespace options {

extern BooleanOptionKey const options; // Show the user-specified options and values
extern BooleanOptionKey const user; // Show the user-specified options and values
extern BooleanOptionKey const all; // Show all the options and values

namespace table {

extern BooleanOptionKey const table;
extern BooleanOptionKey const text; // Show the option definitions table in text format
extern BooleanOptionKey const Wiki; // Show the option definitions table in Wiki format

} // namespace table

extern BooleanOptionKey const exit; // Exit after displaying the options

} // namespace options


/// @brief Lookup functors
typedef  OptionKey  KeyType;
#include <utility/keys/KeyLookup.functors.hh>


} // namespace OptionKeys
} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_keys_OptionKeys_HH
