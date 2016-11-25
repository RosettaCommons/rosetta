// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/qsar/scoring_grid/schema_util.cc
/// @brief  Utility functions for defining XML Schema for scoring grids
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/qsar/scoring_grid/schema_util.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {

/// @brief Used to name the xs:complexType for a scoring grid that is
/// created with the given element name
std::string
complex_type_name_for_scoring_grid( std::string const & element_name )
{
	return "scoring_grid_" + element_name + "_type";
}

/// @brief Define the XML schema definition for a scoring grid that has no
/// subelements but does have a set of attributes (aka options).
void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & scoring_grid_name,
	std::string const & description,
	utility::tag::AttributeList const & attributes
)
{


	//Add weight attribute if not already defined

	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_scoring_grid )
		.element_name( scoring_grid_name )
		.description( description )
		.add_attributes( attributes );
	if ( ! utility::tag::attribute_w_name_in_attribute_list( "weight", attributes ) ) {
		ct_gen.add_attribute( utility::tag::XMLSchemaAttribute::required_attribute( "weight", utility::tag::xsct_real, "XRW TO DO" ) );
	}
	if ( ! utility::tag::attribute_w_name_in_attribute_list( "name", attributes ) ) {
		ct_gen.add_optional_name_attribute();
	}
	ct_gen.write_complex_type_to_schema( xsd );
}

}
}
}

