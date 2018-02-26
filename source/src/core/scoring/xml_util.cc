// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/xml_util.cc
/// @brief Score function utilities useful in RosettaScripts.
/// @author Sarel Fleishman (sarelf@u.washington.edu)
///         Jacob Corn (jecorn@u.washington.edu)
///         Rocco Moretti (rmoretti@u.washington.edu)
///         Eva-Maria Strauch (evas01@uw.edu)

// Unit headers
#include <core/scoring/xml_util.hh>

// Core headers
#include <core/types.hh>
#include <core/scoring/util.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <fstream>

static basic::Tracer TR( "core.scoring.xml_util" );

namespace core {
namespace scoring {

using namespace std;
using namespace core;
using utility::vector1;

/// @details Utility function to find a scorefunction from parser-provided
/// data.
core::scoring::ScoreFunctionOP
parse_score_function(
	utility::tag::TagCOP tag,
	std::string const & option_name,
	basic::datacache::DataMap const & data,
	std::string const & dflt_key/*="commandline"*/ )
{
	std::string const scorefxn_key( tag->getOption<std::string>(option_name, dflt_key) );
	if ( ! data.has( "scorefxns", scorefxn_key ) ) {
		if ( ! core::scoring::check_score_function_sanity( basic::options::option, scorefxn_key ) ) {
			// Incompatible default score functions aren't loaded.
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "ScoreFunction " + scorefxn_key + " was requested with incompatible options.");
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "ScoreFunction " + scorefxn_key + " not found in basic::datacache::DataMap. To add a score function to the data map, define a score function in the <SCOREFXNS/>.'");
		}
	}

	// We don't check "sanity" for the scorefxn_key, as that's user defined, and not necessarily representative of the actual scorefunction.

	return data.get_ptr< ScoreFunction >( "scorefxns", scorefxn_key );
}

/// @details Utility function to find a scorefunction from
/// parser-provided data for the option 'scorefxn'.
core::scoring::ScoreFunctionOP
parse_score_function(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	std::string const & dflt_key/*="commandline"*/ )
{
	return parse_score_function(tag, "scorefxn", data, dflt_key);
}

std::string
get_score_function_name(
	utility::tag::TagCOP tag,
	std::string const & option_name
) {
	return tag->getOption<std::string>(option_name, "commandline");
}

std::string
get_score_function_name(
	utility::tag::TagCOP tag
) {
	return get_score_function_name(tag, "scorefxn");
}

///  Get attributes ( i.e. options ) for movers to build xml schemas
void
attributes_for_get_score_function_name( utility::tag::AttributeList & attributes )
{
	attributes_for_get_score_function_name_w_description( attributes, "" );
}

void
attributes_for_get_score_function_name(
	utility::tag::AttributeList & attributes,
	std::string const & option_name )
{
	attributes_for_get_score_function_name_w_description( attributes, option_name, "" );
}

void
attributes_for_get_score_function_name_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & description )
{
	attributes_for_get_score_function_name_w_description( attributes, "scorefxn", description);
}

void
attributes_for_get_score_function_name_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & option_name,
	std::string const & description )
{
	attributes + utility::tag::XMLSchemaAttribute( option_name, utility::tag::xs_string,
		( description == "" ? "Name of score function to use" : description ) );
}

void
attributes_for_parse_score_function( utility::tag::AttributeList & attributes )
{
	attributes_for_parse_score_function_w_description( attributes, "" );
}

void
attributes_for_parse_score_function( utility::tag::AttributeList & attributes, std::string const & sfxn_option_name )
{
	attributes_for_parse_score_function_w_description( attributes, sfxn_option_name, "" );
}

void
attributes_for_parse_score_function_w_description(
	utility::tag::AttributeList & attributes,
	std::string const & description )
{
	attributes_for_parse_score_function_w_description( attributes, "scorefxn", description );
}

void
attributes_for_parse_score_function_w_description( utility::tag::AttributeList & attributes,
	std::string const & sfxn_option_name,
	std::string const & description )
{
	attributes + utility::tag::XMLSchemaAttribute(
		sfxn_option_name, utility::tag::xs_string ,
		( description == "" ? "Name of score function to use" : description ) );
}

void
attributes_for_parse_score_function_w_description_when_required( utility::tag::AttributeList & attributes,
	std::string const & sfxn_option_name,
	std::string const & description
) {
	attributes + utility::tag::XMLSchemaAttribute::required_attribute(
		sfxn_option_name, utility::tag::xs_string,
		( description == "" ? "Name of score function to use" : description ) );
}



}
}
