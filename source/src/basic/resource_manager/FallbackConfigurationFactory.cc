// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/FallbackConfigurationFactory.cc
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

//unit headers
#include <basic/resource_manager/FallbackConfigurationFactory.hh>

// package headers
#include <basic/resource_manager/FallbackConfiguration.hh>
#include <basic/resource_manager/FallbackConfigurationCreator.hh>

// utility headers
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// Singleton instance and mutex static data members
namespace utility {

using basic::resource_manager::FallbackConfigurationFactory;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< FallbackConfigurationFactory >::singleton_mutex_{};
template <> std::atomic< FallbackConfigurationFactory * > utility::SingletonBase< FallbackConfigurationFactory >::instance_( 0 );
#else
template <> FallbackConfigurationFactory * utility::SingletonBase< FallbackConfigurationFactory >::instance_( 0 );
#endif

}

namespace basic {
namespace resource_manager {

FallbackConfigurationFactory *
FallbackConfigurationFactory::create_singleton_instance()
{
	return new FallbackConfigurationFactory;
}

FallbackConfigurationFactory::FallbackConfigurationFactory() :
	throw_on_double_registration_( false )
{}

FallbackConfigurationOP
FallbackConfigurationFactory::create_fallback_configuration( std::string const & resource_description ) const
{
	FallbackConfigurationCreatorsMap::const_iterator iter = creators_map_.find( resource_description );
	if ( iter == creators_map_.end() ) {
		throw utility::excn::EXCN_Msg_Exception( "No FallbackConfigurationCreator resposible for the FallbackConfiguration named " + resource_description + " was found in the FallbackConfigurationFactory.  Was it correctly registered?" );
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
			throw utility::excn::EXCN_Msg_Exception( errmsg );
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
	FallbackConfigurationCreatorsMap::const_iterator resources( creators_map_.find( desc ));
	return (resources != creators_map_.end());
}

} // namespace resource_manager
} // namespace basic

