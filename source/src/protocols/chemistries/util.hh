// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/chemistries/util.hh
/// @brief  Utility functions for the chemistry classes, primarily, in constructing XML-Schema type definitions
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_chemistries_util_HH
#define INCLUDED_protocols_chemistries_util_HH

// Core headers
#include <core/types.hh>

// Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace chemistries {

/// @brief Used to name the xs:complexType for a chemistry that is
/// created with the "chem_type" tag-name.  Does so by prepending "chemistry_" and
/// appending "Type" to the "chem_type".
std::string
complex_type_name_for_chemistry( std::string const & chem_type );

/// @brief Define the XML schema definition for a Chemistry which
/// possibly contains attributes/options, but no sub-tags
void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & chem_type,
	std::string const & chem_description,
	utility::tag::AttributeList const & attributes
);


/// @brief Define the XML schema definition for a Chemistry that
/// contains subtags and attributes (aka options).
void
xsd_type_definition_w_attributes_and_repeatable_subelements(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & chem_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes,
	utility::tag::XMLSchemaSimpleSubelementList const & subelements
);

}
}

#endif
