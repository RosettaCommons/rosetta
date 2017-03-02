// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/AndJumpSelector.cc
/// @brief  The AndJumpSelector combines logic from multiple JumpSelectors
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/select/jump_selector/AndJumpSelector.hh>
#include <core/select/jump_selector/AndJumpSelectorCreator.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/jump_selector/JumpSelectorFactory.hh>
#include <core/select/jump_selector/util.hh>

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


AndJumpSelector::AndJumpSelector() {}

/// @brief Copy constructor
///
AndJumpSelector::AndJumpSelector( AndJumpSelector const &src) :
	selectors_( src.selectors_ )
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
JumpSelectorOP AndJumpSelector::clone() const { return JumpSelectorOP( new AndJumpSelector(*this) ); }

AndJumpSelector::~AndJumpSelector() {}

AndJumpSelector::AndJumpSelector(JumpSelectorCOP selector1)
{
	add_jump_selector( selector1 );
}

AndJumpSelector::AndJumpSelector( JumpSelectorCOP selector1, JumpSelectorCOP selector2 )
{
	add_jump_selector( selector1 );
	add_jump_selector( selector2 );
}

JumpSubset
AndJumpSelector::apply( core::pose::Pose const & pose ) const
{
	debug_assert( num_selectors() > 0 );

	// make subset neutral for AND operations
	JumpSubset subset( pose.num_jump(), true );
	for ( auto const & rs : selectors_ ) {
		JumpSubset tmp = rs->apply( pose );
		apply_and_to_subset(tmp, subset);
	}
	return subset;
}

void AndJumpSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{
	// grab the comma-separated list of jump selectors that should be ANDed together
	// from the tag, and then grab each of the indicated jump selectors from the datamap.

	std::list< JumpSelectorCOP > local_selectors;
	if ( tag->hasOption("selectors") ) {
		std::string selectors_str;
		try {
			selectors_str = tag->getOption< std::string >( "selectors" );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to access option 'selectors' from AndJumpSelector::parse_my_tag.\n";
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}
		utility::vector1< std::string > selector_names = utility::string_split( selectors_str, ',' );

		for ( std::string const & selector_name : selector_names ) {
			try {
				JumpSelectorCOP selector = datamap.get_ptr< JumpSelector const >( "JumpSelector", selector_name );
				local_selectors.push_back( selector );
			} catch ( utility::excn::EXCN_Msg_Exception & e ) {
				std::stringstream error_msg;
				error_msg << "Failed to find JumpSelector named '" << selector_name << "' from the Datamap from AndJumpSelector::parse_my_tag.\n";
				error_msg << e.msg();
				throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
			}
		}
	} // hasOption selectors

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
		error_msg << "No JumpSelectors given to the AndJumpSelector; AndJumpSelector requires at least one JumpSelector as input\n";
		throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
	}

	for ( auto const & local_selector : local_selectors ) {
		add_jump_selector( local_selector );
	}
}

void AndJumpSelector::add_jump_selector( JumpSelectorCOP selector )
{
	selectors_.push_back(selector);
}

Size AndJumpSelector::num_selectors() const
{
	return selectors_.size();
}

void
AndJumpSelector::apply_and_to_subset(JumpSubset const & newSubset, JumpSubset & existingSubset) const
{
	debug_assert( existingSubset.size() == newSubset.size() );
	for ( Size ii = 1; ii <= existingSubset.size(); ++ii ) {
		existingSubset[ ii ] = existingSubset[ ii ] && newSubset[ ii ];
	}
}

std::string AndJumpSelector::get_name() const {
	return AndJumpSelector::class_name();
}

std::string AndJumpSelector::class_name() {
	return "And";
}

void
AndJumpSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	utility::tag::AttributeList attributes;
	attributes + utility::tag::XMLSchemaAttribute(
		"selectors", utility::tag::xs_string,
		"Comma separated list of selected jumps" );
	xsd_type_definition_w_attributes_and_optional_subselectors(
		xsd, class_name(),
		"The AndJumpSelector combines the output of multiple JumpSelectors using AND "
		"logic, i.e., only jumps selected by ALL contained JumpSelectors will be selected. "
		"JumpSelecters can be pulled in from a DataMap, from subtags (for JumpSelectors "
		"known to the JumpSelectorFactory) or programmatically through %add_jump_selector.",
		attributes );
}


JumpSelectorOP
AndJumpSelectorCreator::create_jump_selector() const {
	return JumpSelectorOP( new AndJumpSelector );
}

std::string
AndJumpSelectorCreator::keyname() const {
	return AndJumpSelector::class_name();
}

void
AndJumpSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	AndJumpSelector::provide_xml_schema( xsd );
}


} //namespace jump_selector
} //namespace select
} //namespace core
