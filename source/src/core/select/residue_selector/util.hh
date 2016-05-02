// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/util.hh
/// @brief  Utility functions for the residue selector classes, primarily, in constructing XML-Schema type definitions
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_select_residue_selector_util_HH
#define INCLUDED_core_select_residue_selector_util_HH

// Core headers
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <string>
#include <list>
#include <algorithm>

namespace core {
namespace select {
namespace residue_selector {

typedef std::list< utility::tag::XMLSchemaAttribute > AttributeList;

/// @brief Used to name the xs:complexType for a residue selector that is
/// created with the "rs_type" tag-name.  Does so by prepending "rs_" and
/// appending "Type" to the "rs_type".  E.g., "rs_AndType" would be the
/// name given to the complexType to describe the format of the
/// AndResidueSelector.
std::string
complex_type_name_for_residue_selector( std::string const & rs_type );

/// @brief Add a type to the XML Schema Definition that might be used
/// in several places.  Available options are:
/// 1) int_cslist -- a comma-separated list of integers
/// This should be moved to src/utility/tag
void
common_simple_types(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & desired_type
);

/// @brief Define the XML schema definition for a ResidueSelector that
/// contains no other ResidueSelectors but may contain some number of
/// attributes (aka options).
void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	AttributeList const & attributes
);

/// @brief Define the XML schema definition for a ResidueSelector that
/// contains a single ResidueSelector in its set of sub-elements
/// (aka sub-tags) and may contain some number of attributes (aka options).
void
xsd_type_definition_w_attributes_and_subselector(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	AttributeList const & attributes
);

/// @brief Define the XML schema definition for a ResidueSelector that
/// contains more than one ResidueSelector in its set of sub-elements
/// (aka sub-tags) and may contain some number of attributes (aka options).
void
xsd_type_definition_w_attributes_and_subselectors(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	core::Size min_occurrence,
	AttributeList const & attributes
);

/// @brief Define the XML schema definition for a ResidueSelector that
/// contains more than one ResidueSelector in its set of sub-elements
/// (aka sub-tags) and may contain some number of attributes (aka options).
void
xsd_type_definition_w_attributes_and_subselectors(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	core::Size min_occurrence,
	core::Size max_occurrence,
	AttributeList const & attributes
);

/// @brief returns a residue selector given a tag and datamap
/// @details Looks for "residue_selector" option in tag
///          If that option isn't found, returns NULL ptr
///          If that option is found, calls get_residue_selector()
core::select::residue_selector::ResidueSelectorCOP
parse_residue_selector( utility::tag::TagCOP tag, basic::datacache::DataMap const & data );

/// @brief returns a residue selector given a selector's name and datamap
/// @details Looks for selector in the datamap
///          Returns a const ptr to the selector
/// @throws utility::excn::EXCN_Msg_Exception if selector is not found in datamap
core::select::residue_selector::ResidueSelectorCOP
get_residue_selector( std::string const & selector_name, basic::datacache::DataMap const & data );

}
}
}

#endif
