// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/datacache/WriteableCacheableDataFactory.cc
/// @brief
/// @author Justin Porter

#include <basic/datacache/WriteableCacheableDataFactory.hh>

#include <basic/datacache/WriteableCacheableDataCreator.hh>
#include <basic/datacache/WriteableCacheableData.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// Singleton instance and mutex static data members
namespace utility {

using basic::datacache::WriteableCacheableDataFactory;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< WriteableCacheableDataFactory >::singleton_mutex_{};
template <> std::atomic< WriteableCacheableDataFactory * > utility::SingletonBase< WriteableCacheableDataFactory >::instance_( 0 );
#else
template <> WriteableCacheableDataFactory * utility::SingletonBase< WriteableCacheableDataFactory >::instance_( 0 );
#endif

}

namespace basic {
namespace datacache {

static thread_local basic::Tracer tr( "basic.datacache.WriteableCacheableDataFactory", t_trace );

WriteableCacheableDataFactory *
WriteableCacheableDataFactory::create_singleton_instance()
{
	return new WriteableCacheableDataFactory;
}

/// @brief add a WriteableCacheableData prototype, using its default type name as the map key
void
WriteableCacheableDataFactory::factory_register( WriteableCacheableDataCreatorOP creator )
{
	if ( !creator ) {
		throw utility::excn::EXCN_NullPointer( "WriteableCacheableDataFactory recieved a null creator pointer." );
	}

	std::string const& data_type = creator->keyname();

	if ( data_creator_map_.find( data_type ) != data_creator_map_.end() ) {
		throw utility::excn::EXCN_BadInput("WriteableCacheableData::factory_register already has a WriteableCachableData creator with name '"
			+ data_type + "'. Conflicting WriteableCacheableData names" );
	}

	data_creator_map_[ data_type ] = creator;
}

/// @brief return new Data instance by key lookup in data_creator_map_
WriteableCacheableDataOP
WriteableCacheableDataFactory::new_data_instance( std::string const & data_type, std::istream &in )
{
	WriteableCacheableDataMap::const_iterator iter( data_creator_map_.find( data_type ) );

	if ( iter == data_creator_map_.end() ) {
		tr.Error << "[ERROR] " << data_type << " was not registered with WritableCacheableDataFactory, and cannot be initialized." << std::endl
			<< "This is probably a result of not being included by protocols/init.WriteableCacheableDataCreators.ihh and "
			<< "protocols/init.WriteableCacheableDataRegistrators.ihh." << std::endl
			<< "Available WriteableCacheableData types: ";
		for (   WriteableCacheableDataMap::const_iterator data_type_it = data_creator_map_.begin();
				data_type_it != data_creator_map_.end(); ++data_type_it ) {
			tr.Error << data_type_it->first << ", ";
		}
		tr.Error << std::endl;

		throw utility::excn::EXCN_BadInput( "Unregistered WriteableCacheableData type '"
			+ data_type + "'." );
	} else if ( ! iter->second ) {
		throw utility::excn::EXCN_BadInput( "WriteableCacheableDataCreatorOP prototype for "
			+ data_type + " was not registered as NULL." );
	}

	return iter->second->create_data( in );
}

} //namespace datacache
} //namespace basic
