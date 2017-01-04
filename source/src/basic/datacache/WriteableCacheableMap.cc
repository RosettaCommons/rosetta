// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/WriteableCacheableMap.hh
/// @brief
/// @author Justin Porter


// unit headers
#include <basic/datacache/WriteableCacheableMap.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {

WriteableCacheableMap::WriteableCacheableMap() : CacheableData() {}

WriteableCacheableMap::WriteableCacheableMap( WriteableCacheableMap const& )= default;

WriteableCacheableMap::~WriteableCacheableMap() = default;

CacheableDataOP WriteableCacheableMap::clone() const {
	return CacheableDataOP( new WriteableCacheableMap(*this) );
}

WriteableCacheableMap::DataMap & WriteableCacheableMap::map() {
	return map_;
}

WriteableCacheableMap::DataMap const & WriteableCacheableMap::map() const {
	return map_;
}

void WriteableCacheableMap::erase( WriteableCacheableDataOP d ) {
	DataMap::const_iterator it = map_.find( d->datatype() );
	if ( it != map_.end() ) {
		map_[ d->datatype() ].erase( d );
	}
}

WriteableCacheableMap::DataMap::const_iterator
WriteableCacheableMap::begin() const {
	return map_.begin();
}

WriteableCacheableMap::DataMap::const_iterator
WriteableCacheableMap::end() const {
	return map_.end();
}

std::set< WriteableCacheableDataOP >& WriteableCacheableMap::operator[]( std::string const& str ){
	return map_[ str ];
}

WriteableCacheableMap::DataMap::const_iterator
WriteableCacheableMap::find( std::string const& str ) const {
	return map_.find( str );
}

bool WriteableCacheableMap::has( WriteableCacheableDataOP data ){
	if ( map_.find( data->datatype() ) != map_.end() ) {
		return ( map_[ data->datatype() ].find( data ) != map_[ data->datatype() ].end() );
	}
	return false;
}

void WriteableCacheableMap::insert( WriteableCacheableDataOP data ){
	map_[ data->datatype() ].insert( data );
}

WriteableCacheableMapOP WriteableCacheableMap::shared_from_this()
{ return utility::pointer::static_pointer_cast<WriteableCacheableMap>( CacheableData::shared_from_this() ); }


} // namespace datacache
} // namespace basic


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
basic::datacache::WriteableCacheableMap::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( map_ ) ); // DataMap
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
basic::datacache::WriteableCacheableMap::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( map_ ); // DataMap
}

SAVE_AND_LOAD_SERIALIZABLE( basic::datacache::WriteableCacheableMap );
CEREAL_REGISTER_TYPE( basic::datacache::WriteableCacheableMap )

CEREAL_REGISTER_DYNAMIC_INIT( basic_datacache_WriteableCacheableMap )
#endif // SERIALIZATION
