// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/filter_schemas.cc
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit header
#include <protocols/filters/filter_schemas.hh>

// Package headers
#include <utility/tag/XMLSchemaGeneration.hh>

namespace protocols {
namespace filters {

std::string
complex_type_name_for_filter( std::string const & filter_name )
{
	return "filter_" + filter_name + "_type";
}

void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & filter_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
)
{
	// If the developer has already added a name attribute, do not add another one.
	utility::tag::AttributeList local_attrs( attributes );
	if ( ! utility::tag::attribute_w_name_in_attribute_list( "name", local_attrs ) ) {
		local_attrs + utility::tag::optional_name_attribute();
	}


	if ( ! utility::tag::attribute_w_name_in_attribute_list( "confidence", local_attrs ) ) {
		local_attrs + utility::tag::XMLSchemaAttribute::attribute_w_default( "confidence", utility::tag::xsct_real, "Probability that the pose will be filtered out if it does not pass this Filter", "1.0" );
	}



	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_filter )
		.element_name( filter_type )
		.description( description )
		.add_attributes( local_attrs )
		.write_complex_type_to_schema( xsd );
}

void
xsd_type_definition_w_attributes_and_repeatable_subelements(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & filter_type,
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

	if ( ! utility::tag::attribute_w_name_in_attribute_list( "confidence", local_attrs ) ) {
		local_attrs + utility::tag::XMLSchemaAttribute::attribute_w_default( "confidence", utility::tag::xsct_real, "Probability that the pose will be filtered out if it does not pass this Filter", "1.0" );
	}

	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_filter )
		.element_name( filter_type )
		.description( description )
		.add_attributes( local_attrs )
		.set_subelements_repeatable( subelements )
		.write_complex_type_to_schema( xsd );
}


}  // namespace filters
}  // namespace protocols
