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

#include <core/select/residue_selector/util.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <sstream>

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
	rs_def.name( rs_type );
	append_name_and_attributes_to_complex_type( attributes, rs_def );

	std::ostringstream oss;
	rs_def.write_definition( 0, oss );
	xsd.add_top_level_element( rs_type, oss.str() );
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
	rs_def.name( rs_type );
	rs_def.type( utility::tag::xsctt_choice );
	rs_def.add_subelement( rsgroup );
	append_name_and_attributes_to_complex_type( attributes, rs_def );

	std::ostringstream oss;
	rs_def.write_definition( 0, oss );
	xsd.add_top_level_element( rs_type, oss.str() );
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
	rs_def.name( rs_type );
	rs_def.type( utility::tag::xsctt_choice );
	rs_def.add_subelement( rsgroup );
	append_name_and_attributes_to_complex_type( attributes, rs_def );

	std::ostringstream oss;
	rs_def.write_definition( 0, oss );
	xsd.add_top_level_element( rs_type, oss.str() );

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
	rs_def.name( rs_type );
	rs_def.type( utility::tag::xsctt_choice );
	rs_def.add_subelement( rsgroup );
	append_name_and_attributes_to_complex_type( attributes, rs_def );

	std::ostringstream oss;
	rs_def.write_definition( 0, oss );
	xsd.add_top_level_element( rs_type, oss.str() );

	// std::ostringstream oss;
	// oss << "<xs:complexType name=\"" << rs_type << "Type mixed=1\">\n";
	// oss << " <xs:group ref=\"ResidueSelector\" minOccurs=\"" << min_occurrence << "\" maxOccurs=\"" << max_occurrence << "\"/>\n";
	// append_name_and_attributes_to_oss( attributes, oss );
	// oss << "</xs:complexType>\n";
	// return oss.str();
}

}
}
}
