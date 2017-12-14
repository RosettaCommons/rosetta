// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/jump_selector/JumpSelectorFactory.cc
/// @brief  Implementation of the class for instantiating arbitrary JumpSelectors
///         from a string --> JumpSelectorCreator map
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/select/jump_selector/JumpSelectorFactory.hh>

// Package headers
#include <core/select/jump_selector/JumpSelector.hh>
#include <core/select/jump_selector/JumpSelectorCreator.hh>
#include <core/select/jump_selector/util.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>

// Boost headers
#include <boost/bind.hpp>

namespace core {
namespace select {
namespace jump_selector {

void
JumpSelectorFactory::factory_register( JumpSelectorCreatorOP creator )
{
	if ( creator_map_.find( creator->keyname() ) != creator_map_.end() ) {
		std::string err_msg = "Factory Name Conflict: Two or more JumpSelectorCreators registered with the name " + creator->keyname();
		if ( throw_on_double_registration_ ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg );
		} else {
			utility_exit_with_message(  err_msg );
		}
	}
	creator_map_[ creator->keyname() ] = creator;
}

bool JumpSelectorFactory::has_type( std::string const & selector_type ) const
{
	return creator_map_.find( selector_type ) != creator_map_.end();
}

/// @brief Get the XML schema for a given jump selector.
/// @details Throws an error if the jump selector is unknown to Rosetta.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
JumpSelectorFactory::provide_xml_schema(
	std::string const &selector_name,
	utility::tag::XMLSchemaDefinition & xsd
) const {
	if ( ! has_type( selector_name ) ) {
		std::string err_msg =  "No JumpSelectorCreator with the name '" + selector_name + "' has been registered with the JumpSelectorFactory";
		throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg );
	}
	auto iter = creator_map_.find( selector_name );
	iter->second->provide_xml_schema( xsd );
}

JumpSelectorOP JumpSelectorFactory::new_jump_selector(
	std::string const & selector_name,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) const
{
	if ( ! has_type( selector_name ) ) {
		std::string err_msg =  "No JumpSelectorCreator with the name '" + selector_name + "' has been registered with the JumpSelectorFactory";
		throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg );
	}
	auto iter = creator_map_.find( selector_name );
	JumpSelectorOP new_selector = iter->second->create_jump_selector();
	new_selector->parse_my_tag( tag, datamap );
	return new_selector;
}

/// @details By convention, the named assigned to each of the complexTypes for JumpSelectors should be
/// what is returned by the function "complex_type_name_for_jump_selector" (declared in
/// core/select/jump_selectors/util.hh) when given the argument returned by that JumpSelector's
/// JumpSelectorCreator's keyname() function. So long as the writing of XML schema for your jump
/// selector is accomplished by the calling the functions in core/select/jump_selectors/util.hh, then
/// this should happen automatically.
void JumpSelectorFactory::define_jump_selector_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	try {
		utility::tag::define_xml_schema_group(
			creator_map_,
			jump_selector_xml_schema_group_name(),
			& complex_type_name_for_jump_selector,
			xsd );
	} catch ( utility::excn::Exception const & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Could not generate an XML Schema for JumpSelectors from JumpSelectorFactory; offending class"
			" must call core::select::jump_selector::complex_type_name_for_jump_selector when defining"
			" its XML Schema\n" + e.msg() );
	}

}

std::map< std::string, JumpSelectorCreatorOP > const &
JumpSelectorFactory::creator_map() const {
	return creator_map_;
}

std::string JumpSelectorFactory::jump_selector_xml_schema_group_name()
{
	return "jump_selector";
}


void JumpSelectorFactory::set_throw_on_double_registration()
{
	throw_on_double_registration_ = true;
}


JumpSelectorFactory::JumpSelectorFactory() :
	throw_on_double_registration_( false )
{}


} //namespace jump_selector
} //namespace select
} //namespace core
