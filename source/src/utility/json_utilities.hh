// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file utility/json_utilities.hh
/// @brief bag of utility functions for nlohmann::json (the one at external/include/json.hpp)
/// @details These utility functions offer "safe access" to json - it checks that the element is present exactly once and returns-by-reference the contents.  This uses Rosetta-type exceptions/exits to do so.
/// @author Steven Lewis (smlewi@gmail.com)

#ifndef INCLUDED_utility_json_utilities_HH
#define INCLUDED_utility_json_utilities_HH

#if defined(__clang__)
#if (__clang_major__ * 10000 + __clang_minor__ * 100 + __clang_patchlevel__) < 30400
#define _NLOHMANN_JSON_DISABLED_
#endif
#elif defined(__GNUC__) && !(defined(__ICC) || defined(__INTEL_COMPILER))
#if (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__) < 40900
#define _NLOHMANN_JSON_DISABLED_
#endif
#endif

#if not defined(_NLOHMANN_JSON_DISABLED_)
#define _NLOHMANN_JSON_ENABLED_
#endif

#ifdef _NLOHMANN_JSON_ENABLED_

// Unit headers
#include <json.hpp> //external/include/json.hpp

// Protocol headers

// Core headers
#include <platform/types.hh>

// Basic/Utility headers
#include <string>

namespace utility {

//NOTE: variables are typed nlohmann::json instead of just json to prevent collision with the namespace name

/// @brief verifies that the element is present exactly once.  returns at all if true; throws utility_exist otherwise
void verify_present_exactly_once_in_json( nlohmann::json const & json, std::string const & element_name );

/// @brief utility function - verifies that a particular boolean is in the json exactly once
void extract_boolean_from_json( nlohmann::json const & json, std::string const & bool_name, bool & return_bool );

/// @brief utility function - verifies that a particular numeric value is in the json exactly once.  Note JSON spec does not discriminate float vs int; we use floating point here, and your ints will be ok.  You will need to cast it back to core::Size yourself after!
void extract_number_from_json( nlohmann::json const & json, std::string const & number_name, platform::Real & return_number );

/// @brief utility function - verifies that a particular string is in the json exactly once, and is nonempty
void extract_nonempty_string_from_json( nlohmann::json const & json, std::string const & string_name, std::string & return_string );

/// @brief utility function - verifies that a particular json array is in the json exactly once, and is nonempty
void extract_nonempty_array_from_json( nlohmann::json const & json, std::string const & array_name, nlohmann::json & return_array );

/// @brief utility function - verifies that a particular json object is in the json exactly once, and is nonempty
void extract_nonempty_object_from_json( nlohmann::json const & json, std::string const & object_name, nlohmann::json & return_object );




/// @brief lookup for element `key` and if it's type is compatible extract it to `value, return `true` on success, no exception is raised if key is not found
bool extract_value_if_present(json const &j, std::string const & key, bool &value);

/// @brief lookup for element `key` and if it's type is compatible extract it to `value, return `true` on success, no exception is raised if key is not found
bool extract_value_if_present(json const &j, std::string const & key, int &value);

/// @brief lookup for element `key` and if it's type is compatible extract it to `value, return `true` on success, no exception is raised if key is not found
bool extract_value_if_present(json const &j, std::string const & key, double &value);

/// @brief lookup for element `key` and if it's type is compatible extract it to `value, return `true` on success, no exception is raised if key is not found
bool extract_value_if_present(json const &j, std::string const & key, std::string &value);

} // utility

#endif // ifdef _NLOHMANN_JSON_ENABLED_

#endif //INCLUDED_utility_json_utilities_HH
