// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/resource_manager/locator_schemas.cc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit header
#include <basic/resource_manager/locator/locator_schemas.hh>

// package headers
#include <basic/resource_manager/ResourceLocatorFactory.hh>

// Package headers
#include <utility/tag/XMLSchemaGeneration.hh>

namespace basic {
namespace resource_manager {
namespace locator {

void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & locator_type,
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
	ct_gen.complex_type_naming_func( & ResourceLocatorFactory::complex_type_name_for_locator )
		.element_name( locator_type )
		.description( description )
		.add_attributes( local_attrs )
		.write_complex_type_to_schema( xsd );
}

void
xsd_type_definition_w_attributes_and_repeatable_subelements(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & locator_type,
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
	ct_gen.complex_type_naming_func( & ResourceLocatorFactory::complex_type_name_for_locator )
		.element_name( locator_type )
		.description( description )
		.add_attributes( local_attrs )
		.set_subelements_repeatable( subelements )
		.write_complex_type_to_schema( xsd );
}

}  // namespace locator
}  // namespace resource_manager
}  // namespace basic
