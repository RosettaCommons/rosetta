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
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>

namespace core {
namespace select {
namespace residue_selector {

#ifdef MULTITHREADED
std::atomic< ResidueSelectorFactory * > ResidueSelectorFactory::instance_( 0 );
#else
ResidueSelectorFactory * ResidueSelectorFactory::instance_( 0 );
#endif

#ifdef MULTI_THREADED

std::mutex ResidueSelectorFactory::singleton_mutex_;

std::mutex & ResidueSelectorFactory::singleton_mutex() { return singleton_mutex_; }

#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
ResidueSelectorFactory * ResidueSelectorFactory::get_instance()
{
	boost::function< ResidueSelectorFactory * () > creator = boost::bind( &ResidueSelectorFactory::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

ResidueSelectorFactory *
ResidueSelectorFactory::create_singleton_instance()
{
	return new ResidueSelectorFactory;
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
	} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
		throw utility::excn::EXCN_Msg_Exception( "Could not generate an XML Schema for ResidueSelectors from ResidueSelectorFactory; offending class"
			" must call core::select::residue_selector::complex_type_name_for_residue_selector when defining"
			" its XML Schema\n" + e.msg() );
	}

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
