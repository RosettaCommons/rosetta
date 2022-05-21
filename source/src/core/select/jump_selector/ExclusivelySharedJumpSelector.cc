// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/ExclusivelySharedJumpSelector.cc
/// @brief  ExclusivelySharedJumpSelector selects the jump that builds ALL and ONLY the residues passed in
/// @details unit tested in EnsureExclusivelySharedJumpMover.cxxtest.hh
/// @author Jack Maguire

// Unit headers
#include <core/select/jump_selector/ExclusivelySharedJumpSelector.hh>
#include <core/select/jump_selector/ExclusivelySharedJumpSelectorCreator.hh>

// Basic Headers
#include <basic/datacache/DataMap.fwd.hh>

// Package headers
#include <core/select/jump_selector/util.hh>
#include <core/select/residue_selector/util.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/types.hh>
#include <core/select/residue_selector/JumpDownstreamSelector.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <set>

namespace core {
namespace select {
namespace jump_selector {

using namespace residue_selector;

ExclusivelySharedJumpSelector::ExclusivelySharedJumpSelector(
	residue_selector::ResidueSelectorCOP selector
):
	JumpSelector(),
	selector_( selector )
{}

JumpSelectorOP ExclusivelySharedJumpSelector::clone() const { return utility::pointer::make_shared< ExclusivelySharedJumpSelector >(*this); }

JumpSubset
ExclusivelySharedJumpSelector::apply( core::pose::Pose const & pose ) const
{
	runtime_assert_string_msg( selector_ != nullptr, "ExclusivelySharedJumpSelector requires a residue selector" );

	utility::vector1< bool > const desired_selection = selector_->apply( pose );

	JumpSubset subset( pose.num_jump(), false );

	for ( core::Size jump_id = 1; jump_id <= pose.num_jump(); ++jump_id ) {
		residue_selector::JumpDownstreamSelector const jd_selector( jump_id );
		utility::vector1< bool > const selected_by_jump_id = jd_selector.apply( pose );

		//perform element-wise comparison for equality
		subset[ jump_id ] = ( desired_selection == selected_by_jump_id );

		//There should not be multiple jumps that meet this condition.
		//Keep looping to make sure, though
	}

	return subset;
}

void
ExclusivelySharedJumpSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) {
	runtime_assert_string_msg( tag->hasOption( "residue_selector" ),
		"residue_selector attribute was not passed to ExclusivelySharedJumpSelector" );

	std::string const selector_str = tag->getOption< std::string >( "residue_selector" );
	selector_ = get_residue_selector( selector_str, datamap );
}

std::string
ExclusivelySharedJumpSelector::get_name() const {
	return ExclusivelySharedJumpSelector::class_name();
}

std::string
ExclusivelySharedJumpSelector::class_name() {
	return "ExclusivelySharedJumpSelector";
}

void
ExclusivelySharedJumpSelector::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {
	using namespace utility::tag;
	AttributeList attributes;

	core::select::residue_selector::attributes_for_parse_residue_selector(
		attributes,
		"residue_selector",
		"This jump selector will select all of the jumps that exclusively build the residues selected by this residue selector"
	);

	xsd_type_definition_w_attributes(
		xsd, class_name(),
		"ExclusivelySharedJumpSelector selects the jump that builds ALL and ONLY the residues passed in. If no jumps meet these criteria, none are selected.",
		attributes );
}


JumpSelectorOP
ExclusivelySharedJumpSelectorCreator::create_jump_selector() const {
	return utility::pointer::make_shared< ExclusivelySharedJumpSelector >();
}

std::string
ExclusivelySharedJumpSelectorCreator::keyname() const {
	return ExclusivelySharedJumpSelector::class_name();
}

void
ExclusivelySharedJumpSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	ExclusivelySharedJumpSelector::provide_xml_schema( xsd );
}

} //namespace jump_selector
} //namespace select
} //namespace core
