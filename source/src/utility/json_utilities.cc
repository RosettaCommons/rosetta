// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file utility/json/json_utilities.cc
/// @brief bag of utility functions for nlohmann::json (the one at external/include/json.hpp)
/// @author Steven Lewis (smlewi@gmail.com)

// Unit headers
#include <utility/json_utilities.hh>


#ifdef _NLOHMANN_JSON_ENABLED_

// Core headers

// Utility headers
#include <utility/exit.hh>

namespace utility {

/// @brief verifies that the element is present exactly once.  returns at all if true; throws utility_exist otherwise
void verify_present_exactly_once_in_json( nlohmann::json const & json, std::string const & element_name ) {
	if ( json.count(element_name) == 0 ) {
		utility_exit_with_message("JSON element " + element_name + " missing");
	} /*else if ( json.count(element_name) != 1 ) {
	utility_exit_with_message("JSON element " + element_name + " present more than once");
	}*/ //this json library overwrites duplicate elements, it can be only 1 or 0
}

/// @brief utility function - verifies that a particular boolean is in the json exactly once
void extract_boolean_from_json( nlohmann::json const & json, std::string const & bool_name, bool & return_bool ){

	verify_present_exactly_once_in_json( json, bool_name );
	if ( !json[bool_name].is_boolean() ) {
		utility_exit_with_message("JSON element " + bool_name + " not of boolean type");
	}

	return_bool = json[bool_name];
	return;
}

/// @brief utility function - verifies that a particular numeric value is in the json exactly once.  Note JSON spec does not discriminate float vs int; we use floating point here, and your ints will be ok.
void extract_number_from_json( nlohmann::json const & json, std::string const & number_name, platform::Real & return_number ){

	verify_present_exactly_once_in_json( json, number_name );
	if ( !json[number_name].is_number() ) {
		utility_exit_with_message("JSON element " + number_name + " not of number type");
	}

	return_number = json[number_name];
	return;
}

/// @brief utility function - verifies that a particular string is in the json exactly once, and is nonempty
void extract_nonempty_string_from_json( nlohmann::json const & json, std::string const & string_name, std::string & return_string ){

	verify_present_exactly_once_in_json(json, string_name);
	if ( !json[string_name].is_string() ) {
		utility_exit_with_message("JSON element " + string_name + " not of string type");
	}

	return_string = json[string_name].get<std::string>();

	if ( return_string.empty() ) {
		utility_exit_with_message(string_name + " empty");
	}
	return;
}

/// @brief utility function - verifies that a particular json array is in the json exactly once, and is nonempty
void extract_nonempty_array_from_json( nlohmann::json const & json, std::string const & array_name, nlohmann::json & return_array ){

	verify_present_exactly_once_in_json(json, array_name);
	if ( !json[array_name].is_array() ) {
		utility_exit_with_message("JSON element " + array_name + " not of array type");
	}

	return_array = json[array_name];

	if ( return_array.empty() ) {
		utility_exit_with_message(array_name + " empty");
	}
	return;
}

/// @brief utility function - verifies that a particular json object is in the json exactly once, and is nonempty
void extract_nonempty_object_from_json( nlohmann::json const & json, std::string const & object_name, nlohmann::json & return_object ){

	verify_present_exactly_once_in_json(json, object_name);
	if ( !json[object_name].is_object() ) {
		utility_exit_with_message("JSON element " + object_name + " not of object type");
	}

	return_object = json[object_name];

	if ( return_object.empty() ) {
		utility_exit_with_message(object_name + " empty");
	}
	return;
}


/// @brief lookup for element `key` and if it's type is compatible extract it to `value, return `true` on success, no exception is raised if key is not found
bool extract_value_if_present(json const &j, std::string const & key, bool &value)
{
	auto it = j.find(key);
	if( it != j.end() ) {
		if( it->is_boolean() ) {
			value = *it;
			return true;
		}
	}
	return false;
}


/// @brief lookup for element `key` and if it's type is compatible extract it to `value, return `true` on success, no exception is raised if key is not found
bool extract_value_if_present(json const &j, std::string const & key, int &value)
{
	auto it = j.find(key);
	if( it != j.end() ) {
		if( it->is_number() ) {
			value = *it;
			return true;
		}
	}
	return false;
}


/// @brief lookup for element `key` and if it's type is compatible extract it to `value, return `true` on success, no exception is raised if key is not found
bool extract_value_if_present(json const &j, std::string const & key, double &value)
{
	auto it = j.find(key);
	if( it != j.end() ) {
		if( it->is_number() ) {
			value = *it;
			return true;
		}
	}
	return false;
}

/// @brief lookup for element `key` and if it's type is compatible extract it to `value, return `true` on success, no exception is raised if key is not found
bool extract_value_if_present(json const &j, std::string const & key, std::string &value)
{
	auto it = j.find(key);
	if( it != j.end() ) {
		if( it->is_string() ) {
			value = *it;
			return true;
		}
	}
	return false;
}


} //utility

#endif // ifdef _NLOHMANN_JSON_ENABLED_
