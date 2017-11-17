// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/ResidueSelectorFactory.cc
/// @brief  Implementation of the class for instantiating arbitrary ResidueSelectors
///         from a string --> ResidueSelectorCreator map
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreator.hh>
#include <core/select/residue_selector/util.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>

// Boost headers
#include <boost/bind.hpp>

namespace core {
namespace select {
namespace residue_selector {

void
ResidueSelectorFactory::factory_register( ResidueSelectorCreatorOP creator )
{
	if ( creator_map_.find( creator->keyname() ) != creator_map_.end() ) {
		std::string err_msg = "Factory Name Conflict: Two or more ResidueSelectorCreators registered with the name " + creator->keyname();
		if ( throw_on_double_registration_ ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg );
		} else {
			utility_exit_with_message(  err_msg );
		}
	}
	creator_map_[ creator->keyname() ] = creator;
}

bool ResidueSelectorFactory::has_type( std::string const & selector_type ) const
{
	return creator_map_.find( selector_type ) != creator_map_.end();
}

/// @brief Get the XML schema for a given residue selector.
/// @details Throws an error if the residue selector is unknown to Rosetta.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
void
ResidueSelectorFactory::provide_xml_schema(
	std::string const &selector_name,
	utility::tag::XMLSchemaDefinition & xsd
) const {
	if ( ! has_type( selector_name ) ) {
		std::string err_msg =  "No ResidueSelectorCreator with the name '" + selector_name + "' has been registered with the ResidueSelectorFactory";
		throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg );
	}
	std::map< std::string, ResidueSelectorCreatorOP >::const_iterator iter = creator_map_.find( selector_name );
	iter->second->provide_xml_schema( xsd );
}

ResidueSelectorOP ResidueSelectorFactory::new_residue_selector(
	std::string const & selector_name,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) const
{
	if ( ! has_type( selector_name ) ) {
		std::string err_msg =  "No ResidueSelectorCreator with the name '" + selector_name + "' has been registered with the ResidueSelectorFactory";
		throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg );
	}
	std::map< std::string, ResidueSelectorCreatorOP >::const_iterator iter = creator_map_.find( selector_name );
	ResidueSelectorOP new_selector = iter->second->create_residue_selector();
	new_selector->parse_my_tag( tag, datamap );
	return new_selector;
}

/// @details By convention, the named assigned to each of the complexTypes for ResidueSelectors should be
/// what is returned by the function "complex_type_name_for_residue_selector" (declared in
/// core/select/residue_selectors/util.hh) when given the argument returned by that ResidueSelector's
/// ResidueSelectorCreator's keyname() function. So long as the writing of XML schema for your residue
/// selector is accomplished by the calling the functions in core/select/residue_selectors/util.hh, then
/// this should happen automatically.
void ResidueSelectorFactory::define_residue_selector_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	try {
		utility::tag::define_xml_schema_group(
			creator_map_,
			residue_selector_xml_schema_group_name(),
			& complex_type_name_for_residue_selector,
			xsd );
	} catch ( utility::excn::Exception const & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "Could not generate an XML Schema for ResidueSelectors from ResidueSelectorFactory; offending class"
			" must call core::select::residue_selector::complex_type_name_for_residue_selector when defining"
			" its XML Schema\n" + e.msg() );
	}

}

std::map< std::string, ResidueSelectorCreatorOP > const &
ResidueSelectorFactory::creator_map() const {
	return creator_map_;
}

std::string ResidueSelectorFactory::residue_selector_xml_schema_group_name()
{
	return "residue_selector";
}


void ResidueSelectorFactory::set_throw_on_double_registration()
{
	throw_on_double_registration_ = true;
}


ResidueSelectorFactory::ResidueSelectorFactory() :
	throw_on_double_registration_( false )
{}


} //namespace residue_selector
} //namespace select
} //namespace core
