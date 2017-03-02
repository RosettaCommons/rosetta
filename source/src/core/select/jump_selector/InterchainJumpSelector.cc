// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/InterchainJumpSelector.hh
/// @brief  The InterchainJumpSelector selects jumps that land on different chains than they took off from
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/select/jump_selector/InterchainJumpSelector.hh>
#include <core/select/jump_selector/InterchainJumpSelectorCreator.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/jump_selector/util.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/selection.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <utility/assert.hh>

namespace core {
namespace select {
namespace jump_selector {

InterchainJumpSelector::InterchainJumpSelector()
{}

InterchainJumpSelector::~InterchainJumpSelector() {}

/// @brief Copy constructor
InterchainJumpSelector::InterchainJumpSelector( InterchainJumpSelector const & )
{}

JumpSelectorOP InterchainJumpSelector::clone() const { return JumpSelectorOP( new InterchainJumpSelector(*this) ); }

JumpSubset
InterchainJumpSelector::apply( core::pose::Pose const & pose ) const
{
	JumpSubset subset( pose.num_jump(), false );

	kinematics::FoldTree const & ft( pose.fold_tree() );
	for ( core::Size ii = 1; ii <= pose.num_jump(); ++ii ) {
		int ii_up = ft.upstream_jump_residue( ii );
		int ii_down = ft.downstream_jump_residue( ii );
		if ( pose.residue( ii_up ).chain() != pose.residue( ii_down ).chain() ) {
			subset[ ii ] = true;
		}
	}
	return subset;
}

void
InterchainJumpSelector::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &
)
{
}

std::string InterchainJumpSelector::get_name() const {
	return InterchainJumpSelector::class_name();
}

std::string InterchainJumpSelector::class_name() {
	return "Interchain";
}

void
InterchainJumpSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;

	xsd_type_definition_w_attributes(
		xsd, class_name(),
		"TO DO",
		attributes );
}


JumpSelectorOP
InterchainJumpSelectorCreator::create_jump_selector() const {
	return JumpSelectorOP( new InterchainJumpSelector );
}

std::string
InterchainJumpSelectorCreator::keyname() const {
	return InterchainJumpSelector::class_name();
}

void
InterchainJumpSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	InterchainJumpSelector::provide_xml_schema( xsd );
}

} //namespace jump_selector
} //namespace select
} //namespace core
