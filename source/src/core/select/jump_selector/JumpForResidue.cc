// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/JumpForResidue.cc
/// @brief  JumpForResidue selects the jump that builds the residues passed in
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Unit headers
#include <core/select/jump_selector/JumpForResidue.hh>
#include <core/select/jump_selector/JumpForResidueCreator.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/jump_selector/util.hh>
#include <core/select/residue_selector/util.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/selection.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <utility/assert.hh>

namespace core {
namespace select {
namespace jump_selector {

using namespace residue_selector;

//Anonymous utilities
namespace {

core::Size
get_jump_for_resid(
	core::kinematics::FoldTree const & ft,
	core::Size const resid
){
	if ( ft.is_root( resid ) ) {
		utility_exit_with_message( "Attempting to find the jump that builds a residue that is not built by a jump." );
	}

	core::kinematics::Edge const & edge = ft.get_residue_edge( resid );
	if ( edge.is_jump() ) {
		return ft.get_jump_that_builds_residue( resid );
		//return edge.label(); //Do this instead if this code ever somehow becomes performance critical
	} else {
		core::Size const parent_resid = edge.start();
		return get_jump_for_resid( ft, parent_resid );
		//Infinite recursion is prevented by the opening assert
	}

	runtime_assert( false ); //dead code
	return 0; //for safety
}

}

JumpSelectorOP JumpForResidue::clone() const { return utility::pointer::make_shared< JumpForResidue >(*this); }

JumpSubset
JumpForResidue::apply( core::pose::Pose const & pose ) const
{
	runtime_assert_string_msg( selector_ != nullptr, "JumpForResidue requires a residue selector" );

	utility::vector1< core::Size > const resids =
		selection_positions( selector_->apply( pose ) );

	std::set< core::Size > jumps;
	for ( core::Size const resid : resids ) {
		jumps.insert( get_jump_for_resid( pose.fold_tree(), resid ) );
	}

	if ( jumps.size() > 1 && ! allow_multiple_results_ ) {
		utility_exit_with_message( "JumpForResidue was called with residues built from multiple jumps." );
	}

	JumpSubset subset( pose.num_jump(), false );
	for ( core::Size const j : jumps ) subset[ j ] = true;
	return subset;
}

void
JumpForResidue::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) {
	runtime_assert_string_msg( tag->hasOption( "residue_selector" ), "residue_selector attribute was not passed to JumpForResidue" );

	std::string const selector_str =
		tag->getOption< std::string >( "residue_selector" );

	selector_ = get_residue_selector( selector_str, datamap );

	if ( tag->hasOption( "allow_multiple_results" ) ) {
		set_allow_multiple_results( tag->getOption< bool >( "allow_multiple_results" ));
	}

}

std::string
JumpForResidue::get_name() const {
	return JumpForResidue::class_name();
}

std::string
JumpForResidue::class_name() {
	return "JumpForResidue";
}

void
JumpForResidue::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {
	using namespace utility::tag;
	AttributeList attributes;

	core::select::residue_selector::attributes_for_parse_residue_selector(
		attributes,
		"residue_selector",
		"This selector will select all of the jumps that build the residues selected by the residue selector"
	);

	XMLSchemaAttribute allow_multiple_results_attr = XMLSchemaAttribute::attribute_w_default( "allow_multiple_results", xsct_rosetta_bool, "If enabled, allow the selector to select more than one jump", "false" );
	attributes.push_back( allow_multiple_results_attr );

	xsd_type_definition_w_attributes(
		xsd, class_name(),
		"JumpForResidue selects the jump that builds the residues passed in.",
		attributes );
}


JumpSelectorOP
JumpForResidueCreator::create_jump_selector() const {
	return utility::pointer::make_shared< JumpForResidue >();
}

std::string
JumpForResidueCreator::keyname() const {
	return JumpForResidue::class_name();
}

void
JumpForResidueCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	JumpForResidue::provide_xml_schema( xsd );
}

} //namespace jump_selector
} //namespace select
} //namespace core
