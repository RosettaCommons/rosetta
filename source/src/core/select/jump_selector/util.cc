// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/util.cc
/// @brief  Utility functions for the jump selector classes, primarily, in constructing XML-Schema type definitions
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/select/jump_selector/util.hh>

// Package headers
#include <core/select/jump_selector/JumpSelectorFactory.hh>
#include <core/select/jump_selector/JumpSelector.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <sstream>

static basic::Tracer TR( "core.select.jump_selector.util" );

namespace core {
namespace select {
namespace jump_selector {

// Attribute::Attribute() : required( false ) {}
//
// Attribute::Attribute(
//  std::string const & name_in,
//  std::string const & type_in,
//  bool required_in
// ) :
//  name( name_in ),
//  type( type_in ),
//  required( required_in )
// {}

std::string
complex_type_name_for_jump_selector( std::string const & js_type )
{
	return "js_" + js_type + "_type";
}

void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & js_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
)
{
	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_jump_selector )
		.element_name( js_type )
		.description( description )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );
}

void
xsd_type_definition_w_attributes_and_optional_subselector(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & js_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
)
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList subelement;
	subelement.add_group_subelement( & JumpSelectorFactory::jump_selector_xml_schema_group_name );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_jump_selector )
		.element_name( js_type )
		.description( description )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.set_subelements_single_appearance_optional( subelement )
		.write_complex_type_to_schema( xsd );
}

void
xsd_type_definition_w_attributes_and_optional_subselectors(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & js_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
)
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList subelement;
	subelement.add_group_subelement( & JumpSelectorFactory::jump_selector_xml_schema_group_name );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_jump_selector )
		.element_name( js_type )
		.description( description )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subelement )
		.write_complex_type_to_schema( xsd );
}

void
xsd_type_definition_w_attributes_and_optional_subselectors(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & js_type,
	std::string const & description,
	core::Size min_occurrence,
	core::Size max_occurrence,
	utility::tag::AttributeList const & attributes
)
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList subelement;
	subelement.add_group_subelement( & JumpSelectorFactory::jump_selector_xml_schema_group_name );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_jump_selector )
		.element_name( js_type )
		.description( description )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subelement, min_occurrence, max_occurrence )
		.write_complex_type_to_schema( xsd );

}

JumpSelectorCOP
parse_jump_selector( utility::tag::TagCOP tag, basic::datacache::DataMap const & data, std::string const & option_name )
{
	std::string const selectorname = tag->getOption< std::string >( option_name, "" );
	if ( selectorname.empty() ) {
		return JumpSelectorCOP();
	}
	return get_jump_selector( selectorname, data );
}

/// @brief Companion function for parse_jump_selector
/// @brief This assumes the default jump selector option name ("jump_selector").
void
attributes_for_parse_jump_selector_default_option_name(
	utility::tag::AttributeList & attlist,
	std::string const & documentation_string
) {
	attributes_for_parse_jump_selector( attlist, "jump_selector", documentation_string);
}

void
attributes_for_parse_jump_selector(
	utility::tag::AttributeList & attlist,
	std::string const & option_name /* = "jump_selector" */,
	std::string const & documentation_string /* = "" */
)
{
	using namespace utility::tag;
	attlist + XMLSchemaAttribute( option_name, xs_string, documentation_string == "" ? "The name of the already defined JumpSelector that will be used by this object" : documentation_string );
}

void
attributes_for_parse_jump_selector_when_required(
	utility::tag::AttributeList & attlist,
	std::string const & option_name /* = "jump_selector"*/,
	std::string const & documentation_string /* = ""*/
)
{
	using namespace utility::tag;
	attlist + XMLSchemaAttribute::required_attribute( option_name, xs_string, documentation_string == "" ? "The name of the already defined JumpSelector that will be used by this object" : documentation_string );
}

JumpSelectorCOP
get_jump_selector( std::string const & selector_name, basic::datacache::DataMap const & data )
{
	core::select::jump_selector::JumpSelectorCOP selector;
	try {
		selector = data.get_ptr< core::select::jump_selector::JumpSelector const >( "JumpSelector", selector_name );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::stringstream error_msg;
		error_msg << "Failed to find JumpSelector named '" << selector_name << "' in the DataMap.\n";
		error_msg << e.msg();
		throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
	}
	debug_assert( selector );
	TR << "Found jump selector " << selector_name << std::endl;
	return selector;
}


}
}
}
