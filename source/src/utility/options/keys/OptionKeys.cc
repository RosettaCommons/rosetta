// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/keys/OptionKeys.cc
/// @brief  utility::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note   Sample/starter OptionKey collection: Copy this file
///         to your project and adapt to the project namespace
///         and add the project OptionKeys


// Unit headers
#include <utility/options/keys/OptionKeys.hh>


namespace utility {
namespace options {
namespace OptionKeys {


/// @brief Help option keys
BooleanOptionKey const help( "help" ); // Show program help and exit


/// @brief Option display option keys
namespace options {

BooleanOptionKey const options( "options" );
BooleanOptionKey const user( "options:user" ); // Show the user-specified options and values
BooleanOptionKey const all( "options:all" ); // Show all the options and values

namespace table {

BooleanOptionKey const table( "options:table" );
BooleanOptionKey const text( "options:table:text" ); // Generate the option definitions table in text format
BooleanOptionKey const Wiki( "options:table:Wiki" ); // Generate the option definitions table in Wiki format

} // namespace table

BooleanOptionKey const exit( "options:exit" ); // Exit after displaying the options

} // namespace options


/// @brief Lookup functors
#include <utility/keys/KeyLookup.functors.ii>


} // namespace OptionKeys
} // namespace options
} // namespace utility
