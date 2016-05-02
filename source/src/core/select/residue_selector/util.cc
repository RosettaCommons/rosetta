// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/util.cc
/// @brief  Utility functions for the residue selector classes, primarily, in constructing XML-Schema type definitions
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit headers
#include <core/select/residue_selector/util.hh>

// Core headers
#include <core/select/residue_selector/ResidueSelector.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <sstream>

static THREAD_LOCAL basic::Tracer TR( "core.select.residue_selector.util" );

namespace core {
namespace select {
namespace residue_selector {

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
complex_type_name_for_residue_selector( std::string const & rs_type )
{
	return "rs_" + rs_type + "Type";
}

void
common_simple_types(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & desired_type
) {
	if ( desired_type == "int_cslist" ) {
		// std::ostringstream oss;
		// oss << "<xs:simpleType name=\"int_cslist\">\n";
		// oss << " <xs:restriction base=\"xs:string\">\n";
		// oss << "  <xs:pattern values=\"[0-9]+(,[0-9]+)*\"/>\n";
		// oss << " </xs:restriction>\n";
		// oss << "</xs:simpleType>\n";
		utility::tag::XMLSchemaRestriction int_cslist;
		int_cslist.name( "int_cslist" );
		int_cslist.base_type( utility::tag::xs_string );
		int_cslist.add_restriction( utility::tag::xsr_pattern, "[0-9]+(,[0-9]+)*" );

		std::ostringstream oss;
		int_cslist.write_definition( 0, oss );
		xsd.add_top_level_element( "int_cslist", oss.str() );
	}
}

// std::string
// xsd_type_from_string( std::string const & instring ) {
//  if ( instring == "string" ) {
//   return "xs:string";
//  } else if ( instring == "real" ) {
//   return "xs:decimal";
//  } else if ( instring == "int" ) {
//   return "xs:integer";
//  } else if ( instring == "bool" ) {
//   return "xs:boolean";
//  } else {
//   return instring;
//  }
// }

// void
// append_name_and_attributes_to_oss(
//  AttributeList const & attributes,
//  std::ostringstream & oss
// )
// {
//  oss << " <xs:attribute name=\"name\" type=\"xs:string\" use=\"required\"/>\n;";
//  for ( AttributeList::const_iterator iter = attributes.begin(); iter != attributes.end(); ++iter ) {
//   oss << " <xs:attribute name=\"" << iter->name << "\" type=\"" << xsd_type_from_string( iter->type ) << "\"";
//   if ( iter->required ) {
//    oss << " use=\"required\"";
//   }
//   oss << "/>\n";
//  }
// }

void
append_name_and_attributes_to_complex_type(
	AttributeList const & attributes,
	utility::tag::XMLSchemaComplexType & type_definition
)
{
	using utility::tag::XMLSchemaAttribute;
	type_definition.add_attribute( XMLSchemaAttribute( "name", utility::tag::xs_string ));
	for ( AttributeList::const_iterator iter = attributes.begin(), iter_end = attributes.end();
			iter != iter_end; ++iter ) {
		type_definition.add_attribute( *iter );
	}
}


void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	AttributeList const & attributes
)
{
	utility::tag::XMLSchemaComplexType rs_def;
	std::string complex_type_name = complex_type_name_for_residue_selector( rs_type );
	rs_def.name( complex_type_name );
	append_name_and_attributes_to_complex_type( attributes, rs_def );

	std::ostringstream oss;
	rs_def.write_definition( 0, oss );
	xsd.add_top_level_element( complex_type_name, oss.str() );
}

void
xsd_type_definition_w_attributes_and_subselector(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	AttributeList const & attributes
)
{
	utility::tag::XMLSchemaElementOP rsgroup( new utility::tag::XMLSchemaElement );
	rsgroup->group_name( "residue_selector" ); // the "residue_selector" group is written out by the ResidueSelectorFactory directly.

	utility::tag::XMLSchemaComplexType rs_def;
	std::string complex_type_name = complex_type_name_for_residue_selector( rs_type );
	rs_def.name( complex_type_name );
	rs_def.type( utility::tag::xsctt_choice );
	rs_def.add_subelement( rsgroup );
	append_name_and_attributes_to_complex_type( attributes, rs_def );

	std::ostringstream oss;
	rs_def.write_definition( 0, oss );
	xsd.add_top_level_element( complex_type_name, oss.str() );
}

void
xsd_type_definition_w_attributes_and_subselectors(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	core::Size min_occurrence,
	AttributeList const & attributes
)
{
	utility::tag::XMLSchemaElementOP rsgroup( new utility::tag::XMLSchemaElement );
	rsgroup->group_name( "residue_selector" );
	rsgroup->min_occurs( min_occurrence );
	rsgroup->max_occurs( utility::tag::xsminmax_unbounded );

	utility::tag::XMLSchemaComplexType rs_def;
	std::string complex_type_name = complex_type_name_for_residue_selector( rs_type );
	rs_def.name( complex_type_name );
	rs_def.type( utility::tag::xsctt_choice );
	rs_def.add_subelement( rsgroup );
	append_name_and_attributes_to_complex_type( attributes, rs_def );

	std::ostringstream oss;
	rs_def.write_definition( 0, oss );
	xsd.add_top_level_element( complex_type_name, oss.str() );

	// std::ostringstream oss;
	// oss << "<xs:complexType name=\"" << rs_type << "Type mixed=1\">\n";
	// oss << " <xs:group ref=\"ResidueSelector\" minOccurs=\"" << min_occurrence << "\" maxOccurs=\"unbounded\"/>\n";
	// append_name_and_attributes_to_oss( attributes, oss );
	// oss << "</xs:complexType>\n";
	// return oss.str();
}

void
xsd_type_definition_w_attributes_and_subselectors(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	core::Size min_occurrence,
	core::Size max_occurrence,
	AttributeList const & attributes
)
{
	utility::tag::XMLSchemaElementOP rsgroup( new utility::tag::XMLSchemaElement );
	rsgroup->group_name( "residue_selector" );
	rsgroup->min_occurs( min_occurrence );
	rsgroup->max_occurs( max_occurrence );

	utility::tag::XMLSchemaComplexType rs_def;
	std::string complex_type_name = complex_type_name_for_residue_selector( rs_type );
	rs_def.name( complex_type_name );
	rs_def.type( utility::tag::xsctt_choice );
	rs_def.add_subelement( rsgroup );
	append_name_and_attributes_to_complex_type( attributes, rs_def );

	std::ostringstream oss;
	rs_def.write_definition( 0, oss );
	xsd.add_top_level_element( complex_type_name, oss.str() );

	// std::ostringstream oss;
	// oss << "<xs:complexType name=\"" << rs_type << "Type mixed=1\">\n";
	// oss << " <xs:group ref=\"ResidueSelector\" minOccurs=\"" << min_occurrence << "\" maxOccurs=\"" << max_occurrence << "\"/>\n";
	// append_name_and_attributes_to_oss( attributes, oss );
	// oss << "</xs:complexType>\n";
	// return oss.str();
}

ResidueSelectorCOP
parse_residue_selector( utility::tag::TagCOP tag, basic::datacache::DataMap const & data )
{
	std::string const selectorname = tag->getOption< std::string >( "residue_selector", "" );
	if ( selectorname.empty() ) {
		return ResidueSelectorCOP();
	}
	return get_residue_selector( selectorname, data );
}

ResidueSelectorCOP
get_residue_selector( std::string const & selector_name, basic::datacache::DataMap const & data )
{
	core::select::residue_selector::ResidueSelectorCOP selector;
	try {
		selector = data.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::stringstream error_msg;
		error_msg << "Failed to find ResidueSelector named '" << selector_name << "' in the DataMap.\n";
		error_msg << e.msg();
		throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
	}
	debug_assert( selector );
	TR << "Found residue selector " << selector_name << std::endl;
	return selector;
}


}
}
}
