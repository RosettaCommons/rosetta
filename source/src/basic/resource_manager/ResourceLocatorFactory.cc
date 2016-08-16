// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/ResourceLocatorFactory.cc
/// @brief
/// @author

//unit headers
#include <basic/resource_manager/ResourceLocatorFactory.hh>

// package headers
#include <basic/resource_manager/ResourceLocator.hh>
#include <basic/resource_manager/ResourceLocatorCreator.hh>

// utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
#include <utility/thread/threadsafe_creation.hh>

// Singleton instance and mutex static data members
namespace utility {

using basic::resource_manager::ResourceLocatorFactory;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< ResourceLocatorFactory >::singleton_mutex_{};
template <> std::atomic< ResourceLocatorFactory * > utility::SingletonBase< ResourceLocatorFactory >::instance_( 0 );
#else
template <> ResourceLocatorFactory * utility::SingletonBase< ResourceLocatorFactory >::instance_( 0 );
#endif

}

namespace basic {
namespace resource_manager {

ResourceLocatorFactory::~ResourceLocatorFactory() {}

ResourceLocatorFactory *
ResourceLocatorFactory::create_singleton_instance()
{
	return new ResourceLocatorFactory;
}

/// @details Create a resource locator from a tags object
/// @input locator_type This is the type of the resource locator, e.g., DatabaseResourceLocator
/// @input locator_tag This is the name of the instance of the the resource locator, e.g., coming from the 'tag' field in the ResourceLocator tag 'stage_1_resfiles'.
/// @input tags this is the tag structure that is parsed by the ResourceLocator to initialize itself.
ResourceLocatorOP
ResourceLocatorFactory::create_resource_locator(
	std::string const & locator_type,
	std::string const & locator_tag,
	utility::tag::TagCOP tags
) const
{
	std::map< std::string, ResourceLocatorCreatorOP >::const_iterator iter = creator_map_.find( locator_type );
	if ( iter == creator_map_.end() ) {
		throw utility::excn::EXCN_Msg_Exception( "No ResourceLocatorCreator resposible for the ResourceLocator named " + locator_type + " was found in the ResourceLocatorFactory.  Was it correctly registered?" );
	}
	ResourceLocatorOP locator = iter->second->create_resource_locator();
	locator->locator_tag( locator_tag );
	if ( tags ) {
		locator->parse_my_tag( tags );
	}
	return locator;
}

void
ResourceLocatorFactory::factory_register( ResourceLocatorCreatorOP creator )
{
	std::string locator_type = creator->locator_type();
	std::map< std::string, ResourceLocatorCreatorOP >::const_iterator iter = creator_map_.find( locator_type );
	if ( iter != creator_map_.end() ) {
		std::string errmsg( "Double registration of a ResourceLocatorCreator in the ResourceLocatorFactory, named " + locator_type + ". Are there two registrators for this options object, or have you chosen a previously assigned name to a new resource option?" );
		if ( throw_on_double_registration_ ) {
			throw utility::excn::EXCN_Msg_Exception( errmsg );
		} else {
			utility_exit_with_message( errmsg );
		}
	}
	creator_map_[ locator_type ] = creator;
}

void
ResourceLocatorFactory::set_throw_on_double_registration() { throw_on_double_registration_ = true; }


/// singleton has a private constructor
ResourceLocatorFactory::ResourceLocatorFactory() : throw_on_double_registration_( false ) {}


} // namespace resource_manager
} // namespace basic

