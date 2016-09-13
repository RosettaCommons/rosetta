// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/ResourceManagerFactory.cc
/// @brief
/// @author

//unit headers
#include <basic/resource_manager/ResourceManagerFactory.hh>

// package headers
#include <basic/resource_manager/ResourceManagerCreator.hh>
#include <basic/resource_manager/ResourceManager.hh>

// utility headers
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace basic {
namespace resource_manager {

ResourceManagerFactory *
ResourceManagerFactory::create_singleton_instance()
{
	return new ResourceManagerFactory;
}

ResourceManagerFactory::ResourceManagerFactory() {}

ResourceManager *
ResourceManagerFactory::create_resource_manager_from_options_system() const
{
	// temporary: always instantiate a JD2 resource manager
	auto iter = creators_map_.find( "JD2ResourceManager" );
	if ( iter == creators_map_.end() ) {
		utility_exit_with_message( "failed to find JD2ResourceManager; must not have been properly registered" );
	}
	return iter->second->create_resource_manager();
}

void
ResourceManagerFactory::factory_register( ResourceManagerCreatorOP creator )
{
	std::string manager_type = creator->manager_type();
	std::map< std::string, ResourceManagerCreatorOP >::const_iterator iter = creators_map_.find( manager_type );
	if ( iter != creators_map_.end() ) {
		throw utility::excn::EXCN_Msg_Exception( "Double registration of a ResourceManagerCreator in the ResourceManagerFactory, named " + manager_type + ". Are there two registrators for this ResourceManager, or have you chosen a previously assigned name to a new ResourceManager?" );
	}
	creators_map_[ manager_type ] = creator;
}

} // namespace resource_manager
} // namespace basic

