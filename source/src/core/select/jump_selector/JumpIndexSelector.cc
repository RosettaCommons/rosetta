// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/JumpIndexSelector.hh
/// @brief  The JumpIndexSelector selects a jump with a given index
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/select/jump_selector/JumpIndexSelector.hh>
#include <core/select/jump_selector/JumpIndexSelectorCreator.hh>

// Package headers
#include <core/select/jump_selector/util.hh>

// Project headers
#include <core/pose/selection.hh>
#include <core/kinematics/FoldTree.hh>

// Utility Headers
#include <utility/assert.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace select {
namespace jump_selector {

static basic::Tracer TR("core.select.jump_selector.JumpIndexSelector");

JumpIndexSelector::JumpIndexSelector() = default;

/// @brief Copy constructor
///
JumpIndexSelector::JumpIndexSelector( JumpIndexSelector const & ) = default;

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
JumpSelectorOP JumpIndexSelector::clone() const { return JumpSelectorOP( new JumpIndexSelector(*this) ); }

JumpIndexSelector::JumpIndexSelector( int jump, bool soft_error ):
	jump_( jump ),
	soft_error_( soft_error )
{}

JumpIndexSelector::~JumpIndexSelector() = default;

JumpSubset
JumpIndexSelector::apply( core::pose::Pose const & pose ) const
{
	JumpSubset subset( pose.fold_tree().num_jump(), false );

	if ( jump_ <= 0 || core::Size( jump_ ) > pose.fold_tree().num_jump() ) {
		TR.Warning << "Pose does not contain a jump number of " << jump_ << std::endl;
		if ( soft_error_ ) {
			return subset;
		} else {
			utility_exit_with_message("Attempted to access a jump that does not exist.");
		}
	}

	debug_assert( jump_ > 0 && Size( jump_ ) <= subset.size() );
	subset[ jump_ ] = true;
	return subset;
}

void
JumpIndexSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
)
{
	try {
		jump( tag->getOption< int >( "jump" ) );
	} catch ( utility::excn::EXCN_Msg_Exception & e ) {
		std::stringstream err_msg;
		err_msg << "Failed to access required option 'jump' from JumpIndexSelector::parse_my_tag.\n";
		err_msg << e.msg();
		throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
	}
}

void
JumpIndexSelector::jump( int jump )
{
	jump_ = jump;
}

int
JumpIndexSelector::jump() const
{
	return jump_;
}

std::string JumpIndexSelector::get_name() const {
	return JumpIndexSelector::class_name();
}

std::string JumpIndexSelector::class_name() {
	return "JumpIndex";
}

void
JumpIndexSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes + XMLSchemaAttribute::required_attribute(
		"jump", xs_integer,
		"The integer given for the \"jump\" argument should refer to a Jump that is present in the Pose." );

	xsd_type_definition_w_attributes(
		xsd, class_name(),
		"The JumpIndexSelector sets the positions corresponding to all of the "
		"jumps that are upstream of the indicated jump to true, and all the "
		"other positions to false. This selector is logically equivalent to a "
		"NotSelector applied to the JumpDownstreamSelector for the same jump.",
		attributes );
}


JumpSelectorOP
JumpIndexSelectorCreator::create_jump_selector() const {
	return JumpSelectorOP( new JumpIndexSelector );
}

std::string
JumpIndexSelectorCreator::keyname() const {
	return JumpIndexSelector::class_name();
}

void
JumpIndexSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	JumpIndexSelector::provide_xml_schema( xsd );
}

} //namespace jump_selector
} //namespace select
} //namespace core

