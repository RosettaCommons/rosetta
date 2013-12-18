// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/ResFilterFactory.cc
/// @brief
/// @author ashworth

#include <core/pack/task/operation/ResFilterFactory.hh>
#include <core/pack/task/operation/ResFilterCreator.hh>

// AUTO-REMOVED #include <core/pack/task/operation/ResFilters.hh>

#include <utility/exit.hh> // runtime_assert, utility_exit_with_message
// AUTO-REMOVED #include <utility/tag/Tag.hh>

#include <core/pack/task/operation/ResFilter.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace core {
namespace pack {
namespace task {
namespace operation {

// special singleton functions
// initialize
ResFilterFactory * ResFilterFactory::instance_( 0 );

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex ResFilterFactory::singleton_mutex_;

std::mutex & ResFilterFactory::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
ResFilterFactory * ResFilterFactory::get_instance()
{
	boost::function< ResFilterFactory * () > creator = boost::bind( &ResFilterFactory::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

ResFilterFactory *
ResFilterFactory::create_singleton_instance()
{
	return new ResFilterFactory;
}

ResFilterFactory::ResFilterFactory() {}

ResFilterFactory::~ResFilterFactory(){}

void
ResFilterFactory::factory_register( ResFilterCreatorOP creator )
{
	add_creator( creator );
}

///@brief add a ResFilter prototype, using its default type name as the map key
void
ResFilterFactory::add_creator( ResFilterCreatorOP creator )
{
	runtime_assert( creator );
	filter_creator_map_[ creator->keyname() ] = creator;
}

bool ResFilterFactory::has_type( std::string const & type ) const
{
	return ( filter_creator_map_.find( type ) != filter_creator_map_.end() );
}

///@brief return new ResFilter by key lookup in filter_creator_map_ (new ResFilter parses Tag if provided)
ResFilterOP
ResFilterFactory::newResFilter(
	std::string const & type,
	TagCOP tag /* = boost::shared_ptr< Tag >() */
) const
{
	ResFilterCreatorMap::const_iterator iter( filter_creator_map_.find( type ) );
	if ( iter != filter_creator_map_.end() ) {
		ResFilterOP filter( iter->second->create_res_filter() );
		// parse tag if tag pointer is pointing to one
		if ( tag.get() != NULL ) filter->parse_tag( tag );
		return filter;
	} else {
		utility_exit_with_message( type + " is not known to the ResFilterFactory. Was its ResFilterCreator class registered at initialization?" );
		return NULL;
	}
}

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core
