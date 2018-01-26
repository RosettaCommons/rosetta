// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/resource_manager/locator/locator_schemas.hh
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_basic_resource_manager_locator_locator_schemas_HH
#define INCLUDED_basic_resource_manager_locator_locator_schemas_HH

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace basic {
namespace resource_manager {
namespace locator {

/// @brief Define the XML schema definition for a ResourceLocator that
/// contains no subtags but may contain any number of
/// attributes (aka options).
void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & locator_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
);

/// @brief Define the XML schema definition for a ResourceLocator that
/// contains subtags and attributes (aka options).
void
xsd_type_definition_w_attributes_and_repeatable_subelements(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & locator_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes,
	utility::tag::XMLSchemaSimpleSubelementList const & subelements
);

}  // namespace locator
}  // namespace resource_manager
}  // namespace basic

#endif  // INCLUDED_basic_resource_manager_locator_locator_schemas_HH
