// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/ResidueSelectorFactory.cc
/// @brief  Implementation of the class for instantiating arbitrary ResidueSelectors
///         from a string --> ResidueSelectorCreator map
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit headers
#include <core/pack/task/residue_selector/ResidueSelectorFactory.hh>

// Package headers
#include <core/pack/task/residue_selector/ResidueSelector.hh>
#include <core/pack/task/residue_selector/ResidueSelectorCreator.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>

namespace core {
namespace pack {
namespace task {
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

void ResidueSelectorFactory::set_throw_on_double_registration()
{
	throw_on_double_registration_ = true;
}


ResidueSelectorFactory::ResidueSelectorFactory() :
	throw_on_double_registration_( false )
{}


} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core
