// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/qsar/scoring_grid/schema_util.hh
/// @brief  Utility functions for defining XML Schema for scoring grids
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_qsar_scoring_grid_schema_util_HH
#define INCLUDED_protocols_qsar_scoring_grid_schema_util_HH

// Core headers
#include <core/types.hh>

// Utility headers
#include <protocols/qsar/scoring_grid/GridSet.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

// C++ headers
#include <string>
#include <list>

namespace protocols {
namespace qsar {
namespace scoring_grid {

/// @brief Get a GridSet from the datamap, using the option_name attribute of the tag.
GridSetCOP
parse_grid_set_from_tag( utility::tag::TagCOP tag, basic::datacache::DataMap const & data, std::string const & option_name = "grid_set" );

/// @brief Append the attributes read by parse_grid_set_from_tag()
void
attributes_for_parse_grid_set_from_tag(utility::tag::AttributeList &attributes, std::string const & description = "", std::string const & option_name = "grid_set");

/// @brief Get a GridSet from the datamap, using the option_name attribute of the tag.
/// If a GridSet hasn't been specified, return a nullptr
GridSetCOP
parse_optional_grid_set_from_tag( utility::tag::TagCOP tag, basic::datacache::DataMap const & data, std::string const & option_name = "grid_set" );

/// @brief Append the attributes read by parse_optional_grid_set_from_tag()
void
attributes_for_parse_optional_grid_set_from_tag(utility::tag::AttributeList &attributes, std::string const & description = "", std::string const & option_name = "grid_set");

/// @brief Used to name the xs:complexType for a scoring grid that is
/// created with the given element name
std::string
complex_type_name_for_scoring_grid( std::string const & element_name );

/// @brief Define the XML schema definition for a scoring grid that has no
/// subelements but does have a set of attributes (aka options).
void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & scoring_grid_name,
	std::string const & description,
	utility::tag::AttributeList const & attributes
);


}
}
}

#endif
