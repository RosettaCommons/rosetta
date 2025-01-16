// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/chemistries/util.cc
/// @brief  Utility functions for the chemistry classes, primarily, in constructing XML-Schema type definitions
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <protocols/chemistries/util.hh>

// Package headers

// Basic headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <sstream>

static basic::Tracer TR( "protocols.chemistries.util" );

namespace protocols {
namespace chemistries {

std::string
complex_type_name_for_chemistry( std::string const & chem_type )
{
	return "chemistry_" + chem_type + "_type";
}

void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & chem_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
)
{
	// If the developer has already added a name attribute, do not add another one.
	utility::tag::AttributeList local_attrs( attributes );
	if ( ! utility::tag::attribute_w_name_in_attribute_list( "name", local_attrs ) ) {
		local_attrs + utility::tag::optional_name_attribute();
	}

	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_chemistry )
		.element_name( chem_type )
		.description( description )
		.add_attributes( local_attrs )
		.write_complex_type_to_schema( xsd );
}

void
xsd_type_definition_w_attributes_and_repeatable_subelements(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & chem_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes,
	utility::tag::XMLSchemaSimpleSubelementList const & subelements
)
{
	// If the developer has already added a name attribute, do not add another one.
	utility::tag::AttributeList local_attrs( attributes );
	if ( ! utility::tag::attribute_w_name_in_attribute_list( "name", local_attrs ) ) {
		local_attrs + utility::tag::optional_name_attribute();
	}

	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_chemistry )
		.element_name( chem_type )
		.description( description )
		.add_attributes( local_attrs )
		.set_subelements_repeatable( subelements )
		.write_complex_type_to_schema( xsd );
}

}
}
