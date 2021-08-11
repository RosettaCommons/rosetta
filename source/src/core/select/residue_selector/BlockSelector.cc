// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/BlockSelector.hh
/// @brief  selectes a specified continuous block of previously selected residues
/// @author frankdt (frankdt@email.unc.edu)

// Unit headers
#include <core/select/residue_selector/BlockSelector.hh>
#include <core/select/residue_selector/BlockSelectorCreator.hh>
#include <core/select/util.hh>

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


static basic::Tracer TR( "protocols.pose_sewing.residue_selectors.BlockSelector" );

namespace core {
namespace select {
namespace residue_selector {

/// @brief Constructor.
///
BlockSelector::BlockSelector():
	core::select::residue_selector::ResidueSelector()
{
}

/// @brief Destructor.
///
BlockSelector::~BlockSelector() {}

/// @brief Copy Constructor.  Usually not necessary unless you need deep copying (e.g. OPs)
//BlockSelector::BlockSelector(BlockSelector const & src):
// core::select::residue_selector::ResidueSelector( src )
//{
//}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
BlockSelector::clone() const {

	return utility::pointer::make_shared< BlockSelector >( *this );
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
BlockSelector::ResidueSubset
BlockSelector::apply(
	core::pose::Pose const & pose
) const {
	ResidueSubset output_set( pose.size(), false );

	core::select::residue_selector::ResidueSubset selection( pose.total_residue(), false );
	if ( selector_ != nullptr ) {
		selection = selector_->apply( pose );
	}

	utility::vector1<core::select::residue_selector::ResidueSubset> block_selections = core::select::residue_selector::identify_ss_blocks_vec( selection );
	TR << "Block size: " << block_selections.size() << std::endl;
	if ( block_selections.size() == 1 ) {
		utility_exit_with_message("BlockSelector must take a residue selector that is chopped into SS blocks.");
	}

	core::Size block_number = block_number_;
	if ( inverse_ ) {
		block_number = (block_selections.size() + 1 - block_number_);
		TR << "Block number " << block_number << std::endl;
	}
	return block_selections[block_number];
}

void
BlockSelector::set_block_number(core::Size block_number){
	block_number_ = block_number;
}

void
BlockSelector::set_inverse(bool inverse){
	inverse_ = inverse;
}
void
BlockSelector::set_selector(core::select::residue_selector::ResidueSelectorCOP selector){
	selector_ = selector;
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
BlockSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	if ( tag->hasOption("block_number") ) {
		block_number_ = tag->getOption<core::Size>("block_number");
	}
	if ( tag->hasOption("inverse") ) {
		inverse_ = tag->getOption<bool>("inverse");
	}
	if ( tag->hasOption( "selector" ) ) {
		selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "selector" ), datamap );
	} else {
		utility_exit_with_message("A residue selector chopped into contiguous regions is required!");
	}
}

std::string BlockSelector::get_name() const
{
	return BlockSelector::class_name();
}

std::string BlockSelector::class_name()
{
	return "BlockSelector";
}

void BlockSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "block_number", xsct_non_negative_integer, "Which residue block to select", "1" )
		+ XMLSchemaAttribute::attribute_w_default(  "inverse", xs_boolean, "If true, block 1 is the c-terminal-most block", "false" );
	core::select::residue_selector::attributes_for_parse_residue_selector( attributes, "selector",  "residue selector of which a continuous block is to be selected" );
	core::select::residue_selector::xsd_type_definition_w_attributes( xsd, class_name(), "BlockSelector selects regions of continuous blocks from an input selector.  Make sure to use some other selector to chop the selector into countable blocks.", attributes );

}

core::select::residue_selector::ResidueSelectorOP
BlockSelectorCreator::create_residue_selector() const {
	return utility::pointer::make_shared< BlockSelector >();
}

std::string
BlockSelectorCreator::keyname() const {
	return BlockSelector::class_name();
}

/// @brief Provide XSD information, allowing automatic evaluation of bad XML.
void
BlockSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	BlockSelector::provide_xml_schema( xsd );
}

} //residue_selectors
} //pose_sewing
} //protocols

