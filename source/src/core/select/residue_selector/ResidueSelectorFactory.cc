// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/select/residue_selector/ResidueSelectorFactory.cc
/// @brief  Implementation of the class for instantiating arbitrary ResidueSelectors
///         from a string --> ResidueSelectorCreator map
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreator.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

namespace core {
namespace select {
namespace residue_selector {

ResidueSelectorFactory * ResidueSelectorFactory::instance_( 0 );

ResidueSelectorFactory * ResidueSelectorFactory::get_instance()
{
	if ( ! instance_ ) {
		instance_ = new ResidueSelectorFactory;
	}
	return instance_;
}

void
ResidueSelectorFactory::factory_register( ResidueSelectorCreatorOP creator )
{
	if ( creator_map_.find( creator->keyname() ) != creator_map_.end() ) {
		std::string err_msg = "Factory Name Conflict: Two or more ResidueSelectorCreators registered with the name " + creator->keyname();
		if ( throw_on_double_registration_ ) {
			throw utility::excn::EXCN_Msg_Exception( err_msg );
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

ResidueSelectorOP ResidueSelectorFactory::new_residue_selector(
	std::string const & selector_name,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) const
{
	if ( ! has_type( selector_name ) ) {
		std::string err_msg =  "No ResidueSelectorCreator with the name '" + selector_name + "' has been registered with the ResidueSelectorFactory";
		throw utility::excn::EXCN_Msg_Exception( err_msg );
	}
	std::map< std::string, ResidueSelectorCreatorOP >::const_iterator iter = creator_map_.find( selector_name );
	ResidueSelectorOP new_selector = iter->second->create_residue_selector();
	new_selector->parse_my_tag( tag, datamap );
	return new_selector;
}

void ResidueSelectorFactory::define_residue_selector_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	using namespace utility::tag;

	XMLSchemaComplexType rs_group_type;
	std::string residue_selector_group_name( "residue_selector" );
	rs_group_type.name( residue_selector_group_name );
	XMLSchemaComplexTypeOP rs_group( new XMLSchemaComplexType );
	rs_group->type( xsctt_choice );
	for ( std::map< std::string, ResidueSelectorCreatorOP >::const_iterator
			iter = creator_map_.begin(), iter_end = creator_map_.end();
			iter != iter_end; ++iter ) {
		XMLSchemaElementOP rs_element( new XMLSchemaElement );
		rs_element->name( iter->first );
		rs_element->type_name( iter->first + "Type" );
		rs_group->add_subelement( rs_element );
	}
	rs_group_type.subtype( rs_group );

	std::ostringstream oss;
	rs_group_type.write_definition( 0, oss );
	xsd.add_top_level_element( residue_selector_group_name, oss.str() );

	// Now iterate across all ResidueSelectors in the map, and have each one write their definition to the XSD
	for ( std::map< std::string, ResidueSelectorCreatorOP >::const_iterator
			iter = creator_map_.begin(), iter_end = creator_map_.end();
			iter != iter_end; ++iter ) {
		iter->second->provide_selector_xsd( xsd );
	}

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
