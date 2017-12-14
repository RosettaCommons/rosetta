// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/OrJumpSelector.cc
/// @brief  The OrJumpSelector combines logic from multiple JumpSelectors
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/select/jump_selector/OrJumpSelector.hh>
#include <core/select/jump_selector/OrJumpSelectorCreator.hh>

// Package headers
#include <core/select/jump_selector/JumpSelectorFactory.hh>
#include <core/select/jump_selector/util.hh>
#include <core/pose/Pose.fwd.hh>

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


OrJumpSelector::OrJumpSelector() = default;

/// @brief Copy constructor
///
OrJumpSelector::OrJumpSelector( OrJumpSelector const &src) :
	selectors_( src.selectors_ )
{}

OrJumpSelector::~OrJumpSelector() = default;

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
JumpSelectorOP OrJumpSelector::clone() const { return JumpSelectorOP( new OrJumpSelector(*this) ); }

OrJumpSelector::OrJumpSelector( JumpSelectorCOP selector1, JumpSelectorCOP selector2 )
{
	add_jump_selector( selector1 );
	add_jump_selector( selector2 );
}

JumpSubset
OrJumpSelector::apply( core::pose::Pose const & pose ) const
{
	debug_assert( num_selectors() > 0 );

	// make subset neutral for OR operations
	JumpSubset subset( pose.num_jump(), false );
	for ( auto const & rs : selectors_ ) {
		JumpSubset tmp = rs->apply( pose );
		apply_or_to_subset(tmp, subset);
	}
	return subset;
}

void OrJumpSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{
	// grab the comma-separated list of residue selectors that should be ORed together
	// from the tag, and then grab each of the indicated residue selectors from the datamap.
	std::list< JumpSelectorCOP > local_selectors;
	if ( tag->hasOption( "selectors" ) ) {
		std::string selectors_str;
		try {
			selectors_str = tag->getOption< std::string >( "selectors" );
		} catch ( utility::excn::Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to access option 'selectors' from OrJumpSelector::parse_my_tag.\n";
			error_msg << e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
		}
		utility::vector1< std::string > selector_names = utility::string_split( selectors_str, ',' );

		for ( std::string const & selector_name : selector_names ) {
			try {
				JumpSelectorCOP selector = datamap.get_ptr< JumpSelector const >( "JumpSelector", selector_name );
				local_selectors.push_back( selector );
			} catch ( utility::excn::Exception & e ) {
				std::stringstream error_msg;
				error_msg << "Failed to find JumpSelector named '" << selector_name << "' from the Datamap from OrJumpSelector::parse_my_tag.\n";
				error_msg << e.msg();
				throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
			}
		}
	}
	// add selectors from tags
	for ( auto const & itag : tag->getTags() ) {
		JumpSelectorCOP rs = JumpSelectorFactory::get_instance()->new_jump_selector(
			itag->getName(),
			itag,
			datamap
		);
		local_selectors.push_back( rs );
	}

	if ( local_selectors.empty() ) { //size() == 0 ) {
		std::stringstream error_msg;
		error_msg << "No JumpSelectors given to the OrJumpSelector; OrJumpSelector requires at least one JumpSelector as input\n";
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
	}

	for ( auto const & local_selector : local_selectors ) {
		add_jump_selector( local_selector );
	}
}

void OrJumpSelector::add_jump_selector( JumpSelectorCOP selector )
{
	selectors_.push_back(selector);
}

Size OrJumpSelector::num_selectors() const
{
	return selectors_.size();
}

void
OrJumpSelector::apply_or_to_subset(JumpSubset const & newSubset, JumpSubset & existingSubset) const
{
	debug_assert( existingSubset.size() == newSubset.size() );
	for ( Size ii = 1; ii <= existingSubset.size(); ++ii ) {
		existingSubset[ ii ] = existingSubset[ ii ] || newSubset[ ii ];
	}
}

std::string OrJumpSelector::get_name() const {
	return OrJumpSelector::class_name();
}

std::string OrJumpSelector::class_name() {
	return "Or";
}

void
OrJumpSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes;
	attributes + utility::tag::XMLSchemaAttribute( "selectors", utility::tag::xs_string , "Residue selectors that have been defined elsewhere in the script" );
	xsd_type_definition_w_attributes_and_optional_subselectors( xsd, class_name(),"Selector that takes the logical or of the provided residue selectors", attributes );
}


JumpSelectorOP
OrJumpSelectorCreator::create_jump_selector() const {
	return JumpSelectorOP( new OrJumpSelector );
}

std::string
OrJumpSelectorCreator::keyname() const {
	return OrJumpSelector::class_name();
}

void
OrJumpSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	OrJumpSelector::provide_xml_schema( xsd );
}

} //namespace jump_selector
} //namespace select
} //namespace core

