// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/ConfluentSelector.hh
/// @brief  selector that selects all the residues between those selected by another selector that are not interspersed with those selected by a third selector
/// @author frankdt (frankdt@email.unc.edu)

// Unit headers
#include <core/select/residue_selector/ConfluentSelector.hh>
#include <core/select/residue_selector/ConfluentSelectorCreator.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions

// Project headers
#include <core/pose/Pose.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/pointer/memory.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "protocols.pose_sewing.residue_selectors.ConfluentSelector" );

namespace core {
namespace select {
namespace residue_selector {

/// @brief Constructor.
///
ConfluentSelector::ConfluentSelector():
	core::select::residue_selector::ResidueSelector()
{
}

/// @brief Destructor.
///
ConfluentSelector::~ConfluentSelector() {}

/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
//ConfluentSelector::ConfluentSelector(ConfluentSelector const & src):
// core::select::residue_selector::ResidueSelector( src )
//{
//}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
ConfluentSelector::clone() const {

	return utility::pointer::make_shared< ConfluentSelector >( *this );
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
ConfluentSelector::ResidueSubset
ConfluentSelector::apply(
	core::pose::Pose const & pose
) const {
	ResidueSubset output_set( pose.size(), false );
	core::select::residue_selector::ResidueSubset terminus_selection( pose.total_residue(), false );
	if ( terminus_selector_ != nullptr ) {
		terminus_selection = terminus_selector_->apply( pose );
	}
	core::select::residue_selector::ResidueSubset breaking_selection( pose.total_residue(), false );
	if ( breaking_selector_ != nullptr ) {
		breaking_selection = breaking_selector_->apply( pose );
	}

	std::pair<core::Size,core::Size> working_pair;
	std::set<std::pair<core::Size,core::Size>> pairs;
	bool open_pair = false;
	for ( core::Size current_residue = 1; current_residue <= pose.size(); ++current_residue ) {
		if ( breaking_selection[current_residue] ) {
			if ( open_pair ) {
				pairs.insert(working_pair);
				open_pair = false;
			}
		} else if ( terminus_selection[current_residue] ) {
			if ( !open_pair ) {
				open_pair = true;
				working_pair.first = current_residue;
			}
			working_pair.second = current_residue;
		}
	}
	if ( open_pair ) {
		pairs.insert(working_pair);
	}
	for ( auto current_pair : pairs ) {
		for ( core::Size current_residue = current_pair.first; current_residue <= current_pair.second; ++current_residue ) {
			output_set[current_residue] = true;
		}
	}

	return output_set;
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
ConfluentSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{
	if ( tag->hasOption( "terminus_selector" ) ) {
		terminus_selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "terminus_selector" ), datamap );
	}
	if ( tag->hasOption( "breaking_selector" ) ) {
		breaking_selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "breaking_selector" ), datamap );
	}
}

void
ConfluentSelector::set_terminus_selector(core::select::residue_selector::ResidueSelectorCOP terminus_selector)
{
	terminus_selector_ = terminus_selector;
}
void
ConfluentSelector::set_breaking_selector(core::select::residue_selector::ResidueSelectorCOP breaking_selector)
{
	breaking_selector_ = breaking_selector;
}

std::string ConfluentSelector::get_name() const
{
	return ConfluentSelector::class_name();
}

std::string ConfluentSelector::class_name()
{
	return "ConfluentSelector";
}

void ConfluentSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	//Syntax Example:
	using namespace utility::tag;
	AttributeList attributes;
	//attributes
	core::select::residue_selector::attributes_for_parse_residue_selector( attributes, "terminus_selector",  "This defines the ends of the selection. All gaps between residues selected by this selector will be filled in." );
	core::select::residue_selector::attributes_for_parse_residue_selector( attributes, "breaking_selector",  "This selector operates inversely to the terminus selector. All gaps in the terminus selector containing at least one residue selected by the breaking selector are retained as gaps." );
	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), "selector that selects all the residues between those selected by another selector that are not interspersed with those selected by a third selector. In other words, it fills in all the gaps between the N- and C-most residues of the terminus selector that do not have a residue between them selected by the breaking selector.", attributes );

}

core::select::residue_selector::ResidueSelectorOP
ConfluentSelectorCreator::create_residue_selector() const {
	return utility::pointer::make_shared< ConfluentSelector >();
}

std::string
ConfluentSelectorCreator::keyname() const {
	return ConfluentSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
void
ConfluentSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	ConfluentSelector::provide_xml_schema( xsd );
}

} //residue_selectors
} //pose_sewing
} //protocols

