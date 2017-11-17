// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/FallbackConfigurationFactory.cc
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

//unit headers
#include <basic/resource_manager/FallbackConfigurationFactory.hh>

// package headers
#include <basic/resource_manager/FallbackConfigurationCreator.hh>

// utility headers
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/thread/threadsafe_creation.hh>

namespace basic {
namespace resource_manager {

FallbackConfigurationFactory::FallbackConfigurationFactory() :
	throw_on_double_registration_( false )
{}

FallbackConfigurationOP
FallbackConfigurationFactory::create_fallback_configuration( std::string const & resource_description ) const
{
	auto iter = creators_map_.find( resource_description );
	if ( iter == creators_map_.end() ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "No FallbackConfigurationCreator resposible for the FallbackConfiguration named " + resource_description + " was found in the FallbackConfigurationFactory.  Was it correctly registered?" );
	}
	return iter->second->create_fallback_configuration();
}

void
FallbackConfigurationFactory::factory_register( FallbackConfigurationCreatorOP creator )
{
	std::string resource_description = creator->resource_description();
	FallbackConfigurationCreatorsMap::const_iterator iter = creators_map_.find( resource_description );
	if ( iter != creators_map_.end() ) {
		std::string errmsg( "Double registration of a FallbackConfigurationCreator in the FallbackConfigurationFactory, named " + resource_description + ". Are there two registrators for this FallbackConfigurationCreator, or have you chosen a previously assigned name to a new FallbackConfigurationCreator?" );

		if ( throw_on_double_registration_ ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  errmsg );
		} else {
			utility_exit_with_message( errmsg );
		}
	}
	creators_map_[ resource_description ] = creator;
}

void
FallbackConfigurationFactory::set_throw_on_double_registration() { throw_on_double_registration_ = true; }

bool
FallbackConfigurationFactory::has_fallback_for_resource( std::string const & desc ) const
{
	auto resources( creators_map_.find( desc ));
	return (resources != creators_map_.end());
}

} // namespace resource_manager
} // namespace basic

