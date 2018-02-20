// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

/// @brief Used to name the xs:complexType for a residue selector that is
/// created with the "rs_type" tag-name.  Does so by prepending "rs_" and
/// appending "Type" to the "rs_type".  E.g., "rs_AndType" would be the
/// name given to the complexType to describe the format of the
/// AndResidueSelector.
std::string
complex_type_name_for_residue_selector( std::string const & rs_type );

/// @brief Define the XML schema definition for a ResidueSelector that
/// contains no other ResidueSelectors but may contain some number of
/// attributes (aka options).
void
xsd_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	std::string const & rs_description,
	utility::tag::AttributeList const & attributes
);

/// @brief Define the XML schema definition for a ResidueSelector that
/// contains a single ResidueSelector in its set of sub-elements
/// (aka sub-tags) and may contain some number of attributes (aka options).
void
xsd_type_definition_w_attributes_and_optional_subselector(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	std::string const & rs_description,
	utility::tag::AttributeList const & attributes
);

/// @brief Define the XML schema definition for a ResidueSelector that
/// contains more than one ResidueSelector in its set of sub-elements
/// (aka sub-tags) and may contain some number of attributes (aka options).
void
xsd_type_definition_w_attributes_and_optional_subselectors(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	std::string const & rs_description,
	utility::tag::AttributeList const & attributes
);

/// @brief Define the XML schema definition for a ResidueSelector that
/// contains more than one ResidueSelector in its set of sub-elements
/// (aka sub-tags) and may contain some number of attributes (aka options).
void
xsd_type_definition_w_attributes_and_optional_subselectors(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	std::string const & rs_description,
	core::Size min_occurrence,
	core::Size max_occurrence,
	utility::tag::AttributeList const & attributes
);

/// @brief returns a residue selector given a tag and datamap
/// @details Looks for "residue_selector" (or whatever option_name is set to)
///  option in tag.
///          If that option isn't found, returns NULL ptr
///          If that option is found, calls get_residue_selector()
ResidueSelectorCOP
parse_residue_selector(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap const & data,
	std::string const &option_name="residue_selector"
);

/// @brief Companion function for parse_residue_selector
/// @brief This assumes the default residue selector option name ("residue_selector").
void
attributes_for_parse_residue_selector_default_option_name(
	utility::tag::AttributeList & attlist,
	std::string const & documentation_string = ""
);

/// @brief Companion function for parse_residue_selector
void
attributes_for_parse_residue_selector(
	utility::tag::AttributeList & attlist,
	std::string const & option_name = "residue_selector",
	std::string const & documentation_string = ""
);

/// @brief Companion function for parse_residue_selector to be used when it is unacceptible
/// for the parse_residue_selector function to return a null pointer
void
attributes_for_parse_residue_selector_when_required(
	utility::tag::AttributeList & attlist,
	std::string const & option_name = "residue_selector",
	std::string const & documentation_string = ""
);

/// @brief returns a residue selector given a selector's name and datamap
/// @details Looks for selector in the datamap
///          Returns a const ptr to the selector
/// @throws utility::excn::EXCN_Msg_Exception if selector is not found in datamap
ResidueSelectorCOP
get_residue_selector( std::string const & selector_name, basic::datacache::DataMap const & data );

/// @brief returns a residue selector embeded in a tag.
/// for example, the Chain selector in the example below:
///
/// Example:
///<RESIDUE_SELECTORS>
///    <Neighborhood name="chAB_neighbors">
///        <Chain chains="A,B">
///    </Neighborhood>
///</RESIDUE_SELECTORS>
///
/// @thorws utility::excn::EXCN_Msg_Exception if selector is not embedded ( ! tag->size() > 1)
///
ResidueSelectorOP
get_embedded_residue_selector(  utility::tag::TagCOP tag, basic::datacache::DataMap & datamap );

/// @brief returns residue selectors embeded in a tag.
/// for example, the Chain selector in the example below:
///
/// Example:
///<RESIDUE_SELECTORS>
///    <Neighborhood name="chAB_neighbors">
///        <Chain chains="A,B">
///        <Glycan >
///    </Neighborhood>
///</RESIDUE_SELECTORS>
///
/// @thorws utility::excn::EXCN_Msg_Exception if no embedded selectors (tag->size() <= 1)
///
utility::vector1< ResidueSelectorOP >
get_embedded_residue_selectors( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap );

/// @brief If sele1 is already an OrResidueSelector, add sele2 to it.
/// If not, return a new OrResidueSelector which combines the two.
/// @details If either of the selectors are a nullptr, just return the other.
ResidueSelectorOP
OR_combine( ResidueSelectorOP sele1, ResidueSelectorOP sele2 );

/// @brief If sele1 is already an AndResidueSelector, add sele2 to it.
/// If not, return a new AndResidueSelector which combines the two.
/// @details If either of the selectors are a nullptr, just return the other.
ResidueSelectorOP
AND_combine( ResidueSelectorOP sele1, ResidueSelectorOP sele2 );

/// @brief Returns True if all the positions in the ResidueSubset are False
bool all_false_selection( ResidueSubset const & selection );
/// @brief Returns True if all the positions in the ResidueSubset are True
bool all_true_selection( ResidueSubset const & selection );
/// @brief Returns True if at least one position in the ResidueSubset is False
bool has_any_false_selection( ResidueSubset const & selection );
/// @brief Returns True if at least one position in the ResidueSubset is True
bool has_any_true_selection( ResidueSubset const & selection );
/// @brief Returns the number of selected residues in the ResidueSubset
core::Size count_selected( ResidueSubset const & selection );
/// @brief Returns the Rosetta Numbering corresponding to the selected residues
utility::vector1< core::Size > selection_positions( ResidueSubset const & selection );
/// @brief Evaluate if two ResidueSubsets are equal
bool are_selections_equal( ResidueSubset const & selection1, ResidueSubset const & selection2 );
/// @brief Returns a string representing the ResidueSubset (- for non selected, * for selected)
std::string represent_residue_selector( ResidueSubset const & selection, std::string const & is_true="*", std::string const & is_false="-" );

/// @brief Read a string from the input tag with the indicated attribute name and interpret it
/// as boolean logical operations on a set of existing ResidueSelectors that are already stored
/// in the DataMap.
ResidueSelectorOP
parse_residue_selector_logic_string(
	basic::datacache::DataMap const & dm,
	utility::tag::TagCOP const & tag,
	std::string const & selector_logic_attribute_name = "selector_logic"
);

/// @brief Append an (optional) attribute that will be used in the call to parse_residue_selector_logic_string
void
attributes_for_parse_residue_selector_logic_string(
	utility::tag::XMLSchemaDefinition & xsd,
	utility::tag::AttributeList & attributes,
	std::string const & selector_logic_attribute_name = "selector_logic"
);

}
}
}

#endif
