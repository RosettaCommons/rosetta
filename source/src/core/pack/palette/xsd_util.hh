// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/palette/xsd_util.cc
/// @brief XSD utilities for the PackerPalette class and its sub-classes.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).

// Unit header

// Package headers
#include <utility/tag/XMLSchemaGeneration.hh>

namespace core {
namespace pack {
namespace palette {

/// @brief Construct a type name for a PackerPalette.
std::string complex_type_name_for_packer_palette( std::string const & packer_palette_name );

/// @brief Add a type definition for a PackerPalette.
void xsd_type_definition_w_attributes( utility::tag::XMLSchemaDefinition & xsd, std::string const & packer_palette_type, std::string const & description, utility::tag::AttributeList const & attributes );

/// @brief Add a type definition for a PackerPalette, with sub-elements.
void xsd_type_definition_w_attributes_and_repeatable_subelements( utility::tag::XMLSchemaDefinition & xsd, std::string const & packer_palette_type, std::string const & description, utility::tag::AttributeList const & attributes, utility::tag::XMLSchemaSimpleSubelementList const & subelements );

}  // namespace palette
}  // namespace pack
}  // namespace core
