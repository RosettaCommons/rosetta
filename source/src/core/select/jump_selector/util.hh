// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/util.hh
/// @brief  Utility functions for the jump selector classes, primarily, in constructing XML-Schema type definitions
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_select_jump_selector_util_HH
#define INCLUDED_core_select_jump_selector_util_HH

// Core headers
#include <core/select/jump_selector/JumpSelector.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>
#include <list>
#include <algorithm>

namespace core {
namespace select {
namespace jump_selector {

/// @brief Used to name the xs:complexType for a jump selector that is
/// created with the "js_type" tag-name.  Does so by prepending "js_" and
/// appending "Type" to the "js_type".  E.g., "js_AndType" would be the
/// name given to the complexType to describe the format of the
/// AndJumpSelector.
std::string
complex_type_name_for_jump_selector( std::string const & js_type );

/// @brief Define the XML schema definition for a JumpSelector that
/// contains no other JumpSelectors but may contain some number of
/// attributes (aka options).
void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & js_type,
	std::string const & js_description,
	utility::tag::AttributeList const & attributes
);

/// @brief Define the XML schema definition for a JumpSelector that
/// contains a single JumpSelector in its set of sub-elements
/// (aka sub-tags) and may contain some number of attributes (aka options).
void
xsd_type_definition_w_attributes_and_optional_subselector(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & js_type,
	std::string const & js_description,
	utility::tag::AttributeList const & attributes
);

/// @brief Define the XML schema definition for a JumpSelector that
/// contains more than one JumpSelector in its set of sub-elements
/// (aka sub-tags) and may contain some number of attributes (aka options).
void
xsd_type_definition_w_attributes_and_optional_subselectors(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & js_type,
	std::string const & js_description,
	utility::tag::AttributeList const & attributes
);

/// @brief Define the XML schema definition for a JumpSelector that
/// contains more than one JumpSelector in its set of sub-elements
/// (aka sub-tags) and may contain some number of attributes (aka options).
void
xsd_type_definition_w_attributes_and_optional_subselectors(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & js_type,
	std::string const & js_description,
	core::Size min_occurrence,
	core::Size max_occurrence,
	utility::tag::AttributeList const & attributes
);

/// @brief returns a jump selector given a tag and datamap
/// @details Looks for "jump_selector" (or whatever option_name is set to)
///  option in tag.
///          If that option isn't found, returns NULL ptr
///          If that option is found, calls get_jump_selector()
JumpSelectorCOP
parse_jump_selector(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	std::string const &option_name="jump_selector"
);

/// @brief Companion function for parse_jump_selector
/// @brief This assumes the default jump selector option name ("jump_selector").
void
attributes_for_parse_jump_selector_default_option_name(
	utility::tag::AttributeList & attlist,
	std::string const & documentation_string = ""
);

/// @brief Companion function for parse_jump_selector
void
attributes_for_parse_jump_selector(
	utility::tag::AttributeList & attlist,
	std::string const & option_name = "jump_selector",
	std::string const & documentation_string = ""
);

/// @brief Companion function for parse_jump_selector to be used when it is unacceptible
/// for the parse_jump_selector function to return a null pointer
void
attributes_for_parse_jump_selector_when_required(
	utility::tag::AttributeList & attlist,
	std::string const & option_name = "jump_selector",
	std::string const & documentation_string = ""
);

/// @brief returns a jump selector given a selector's name and datamap
/// @details Looks for selector in the datamap
///          Returns a const ptr to the selector
/// @throws utility::excn::EXCN_Msg_Exception if selector is not found in datamap
JumpSelectorCOP
get_jump_selector( std::string const & selector_name, basic::datacache::DataMap const & data );



}
}
}

#endif
