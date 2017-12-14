// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/NotJumpSelector.cc
/// @brief  The NotJumpSelector negates the logic of its loaded JumpSelector
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit headers
#include <core/select/jump_selector/NotJumpSelector.hh>
#include <core/select/jump_selector/NotJumpSelectorCreator.hh>

// Package headers
#include <core/select/jump_selector/util.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/jump_selector/JumpSelectorFactory.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <utility/assert.hh>

namespace core {
namespace select {
namespace jump_selector {


NotJumpSelector::NotJumpSelector() = default;
NotJumpSelector::~NotJumpSelector() = default;

/// @brief Copy constructor
///
NotJumpSelector::NotJumpSelector( NotJumpSelector const &src) :
	selector_( src.selector_ )
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
JumpSelectorOP NotJumpSelector::clone() const { return JumpSelectorOP( new NotJumpSelector(*this) ); }

NotJumpSelector::NotJumpSelector( JumpSelectorCOP selector )
{
	set_jump_selector( selector );
}

JumpSubset
NotJumpSelector::apply( core::pose::Pose const & pose ) const
{
	debug_assert( selector_ );

	JumpSubset subset = selector_->apply( pose );
	subset.flip();
	return subset;
}

void NotJumpSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{
	if ( tag->hasOption( "selector" ) ) { // fetch selector from datamap

		if ( tag->size() > 1 ) { // has subtags
			throw CREATE_EXCEPTION(utility::excn::Exception,  "NotJumpSelector can negate ONE JumpSelector! Either specify 'selector' option or provide subtags but not BOTH\n" );
		}
		// grab the JumpSelector to be negated from the selector option
		// and then grab each of the indicated jump selectors from the datamap.
		std::string selector_str;
		try {
			selector_str = tag->getOption< std::string >( "selector" );
		} catch ( utility::excn::Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to access required option 'selector' from NotJumpSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
		}

		try {
			JumpSelectorCOP selector = datamap.get_ptr< JumpSelector const >( "JumpSelector", selector_str );
			set_jump_selector(selector);
		} catch ( utility::excn::Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find JumpSelector named '" << selector_str << "' from the Datamap from NotJumpSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
		}
	} else if ( tag->size() > 1 ) { // attempt reading subtag
		utility::vector0< utility::tag::TagCOP > const & tags = tag->getTags();
		if ( tags.size() > 1 ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  "NotJumpSelector takes exactly ONE JumpSelector! Multiple selectors were specified.\n" );
		}
		JumpSelectorCOP rs = JumpSelectorFactory::get_instance()->new_jump_selector(
			tags.front()->getName(),
			tags.front(),
			datamap
		);
		set_jump_selector( rs );
	}

	if ( !selector_ ) {
		std::stringstream error_msg;
		error_msg << "No JumpSelector given to the NotJumpSelector; NotJumpSelector requires a JumpSelector as input\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
	}
}

void NotJumpSelector::set_jump_selector( JumpSelectorCOP selector )
{
	selector_ = selector;
}

std::string NotJumpSelector::get_name() const {
	return NotJumpSelector::class_name();
}

std::string NotJumpSelector::class_name() {
	return "Not";
}

void
NotJumpSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes + XMLSchemaAttribute( "selector", xs_string , "XRW TO DO" );
	xsd_type_definition_w_attributes_and_optional_subselector( xsd, class_name(),"XRW TO DO", attributes );
}


JumpSelectorOP
NotJumpSelectorCreator::create_jump_selector() const {
	return JumpSelectorOP( new NotJumpSelector );
}

std::string
NotJumpSelectorCreator::keyname() const {
	return NotJumpSelector::class_name();
}

void
NotJumpSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	NotJumpSelector::provide_xml_schema( xsd );
}


} //namespace jump_selector
} //namespace select
} //namespace core

