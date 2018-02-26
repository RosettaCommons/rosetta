// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/xml_util.hh
/// @brief Score function utilities useful in RosettaScripts.
/// @author Sarel Fleishman (sarelf@u.washington.edu)
///         Jacob Corn (jecorn@u.washington.edu)
///         Rocco Moretti (rmoretti@u.washington.edu)
///         Eva-Maria Strauch (evas01@uw.edu)

#ifndef INCLUDED_core_scoring_rosetta_scripts_HH
#define INCLUDED_core_scoring_rosetta_scripts_HH

// Unit headers

// Project Headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// C++ headers
#include <string>

namespace core {
namespace scoring {

/// @brief Look up the score function defined in the <SCOREFXNS/>
/// through the given option. Defaults to 'commandline'.
core::scoring::ScoreFunctionOP
parse_score_function(
	utility::tag::TagCOP tag,
	std::string const & option_name,
	basic::datacache::DataMap const & data,
	std::string const & dflt_key="commandline" );

/// @brief Look up the score function defined in the <SCOREFXNS/>
///through the option 'scorefxn='. Defaults to 'commandline'.
core::scoring::ScoreFunctionOP
parse_score_function(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	std::string const & dflt_key="commandline" );

/// @brief Look up the name of assigned score function to the given
///option. Use this to prevent hard coding default score functions into
///protocols.
std::string
get_score_function_name(
	utility::tag::TagCOP tag,
	std::string const & option_name);

/// @brief Look up the name of assigned score function to the 'scorefxn='
///option. Use this to prevent hard coding default score functions into
///protocols.
std::string
get_score_function_name(
	utility::tag::TagCOP tag);


///////////////////// Attributes ///////////////////////////

/// @brief Appends the attributes read by get_score_function_name
void
attributes_for_get_score_function_name(
	utility::tag::AttributeList & attributes );

/// @brief Appends the attributes read by get_score_function_name w/ name argument
void
attributes_for_get_score_function_name(
	utility::tag::AttributeList & attributes,
	std::string const & option_name);

/// @brief Appends the attributes read by get_score_function_name
void
attributes_for_get_score_function_name_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & description );

/// @brief Appends the attributes read by get_score_function_name w/ name argument
void
attributes_for_get_score_function_name_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & option_name,
	std::string const & description );

/// @brief Appends the attributes read by parse_score_function
void
attributes_for_parse_score_function( utility::tag::AttributeList & attributes);

/// @brief Appends the attributes read by parse_score_function w/ name argument
void
attributes_for_parse_score_function( utility::tag::AttributeList & attributes,
	std::string const & sfxn_option_name );

/// @brief Appends the attributes read by parse_score_function with description
void
attributes_for_parse_score_function_w_description( utility::tag::AttributeList & attributes,
	std::string const & description );

/// @brief Appends the attributes read by parse_score_function w/ name argument and description
void
attributes_for_parse_score_function_w_description( utility::tag::AttributeList & attributes,
	std::string const & sfxn_option_name,
	std::string const & description );

/// @brief Appends the attributes read by parse_score_function w/ name argument and description.
/// @details This version appends the attributes as required attributes.
/// @author Vikram K. Mulligan.
void
attributes_for_parse_score_function_w_description_when_required( utility::tag::AttributeList & attributes,
	std::string const & sfxn_option_name,
	std::string const & description = ""
);

}
}


#endif
