// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/ResourceLoaderFactory.cc
/// @brief
/// @author

//unit headers
#include <basic/resource_manager/ResourceLoaderFactory.hh>

// package headers
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceLoaderCreator.hh>

// utility headers
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>

#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace basic {
namespace resource_manager {

ResourceLoaderFactory *
ResourceLoaderFactory::create_singleton_instance()
{
	return new ResourceLoaderFactory;
}

ResourceLoaderOP
ResourceLoaderFactory::create_resource_loader(
	std::string const & loader_type
) const
{
	auto iter = creator_map_.find( loader_type );
	if ( iter == creator_map_.end() ) {
		throw utility::excn::EXCN_Msg_Exception( "No ResourceLoaderCreator resposible for the ResourceLoader named " + loader_type + " was found in the ResourceLoaderFactory.  Was it correctly registered?" );
	}
	return iter->second->create_resource_loader();
}


bool
ResourceLoaderFactory::has_resource_loader(
	std::string const & loader_type
) const
{
	return creator_map_.find( loader_type ) != creator_map_.end();
}

std::list< std::string >
ResourceLoaderFactory::available_resource_loaders() const
{
	std::list< std::string > loader_types;
	for ( auto const & iter : creator_map_ ) {
		loader_types.push_back( iter.first );
	}
	return loader_types;
}

void
ResourceLoaderFactory::factory_register( ResourceLoaderCreatorOP creator )
{
	std::string loader_type = creator->loader_type();
	std::map< std::string, ResourceLoaderCreatorOP >::const_iterator iter = creator_map_.find( loader_type );
	if ( iter != creator_map_.end() ) {
		std::string errmsg("Double registration of a ResourceLoaderCreator in the ResourceLoaderFactory, named " + loader_type + ". Are there two registrators for this ResourceLoader object, or have you chosen a previously assigned name to a new resource option?" );
		if ( throw_on_double_registration_ ) {
			throw utility::excn::EXCN_Msg_Exception( errmsg );
		} else {
			utility_exit_with_message( errmsg );
		}
	}
	creator_map_[ loader_type ] = creator;
}

void
ResourceLoaderFactory::set_throw_on_double_registration() { throw_on_double_registration_ = true; }

/// singleton has a private constructor
ResourceLoaderFactory::ResourceLoaderFactory() : throw_on_double_registration_( false ) {}


} // namespace resource_manager
} // namespace basic

