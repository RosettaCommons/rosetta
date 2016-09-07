// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Sarel Fleishman

// Unit Headers
#include <basic/datacache/ConstDataMap.hh>

// Package headers

// Project headers
// ObjexxFCL Headers

// C++ Headers
// AUTO-REMOVED #include <string>
#include <map>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <basic/Tracer.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/pointer/ReferenceCount.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {

static THREAD_LOCAL basic::Tracer TR( "basic.datacache.ConstDataMap" );

ConstDataMap::ConstDataMap() {}
ConstDataMap::ConstDataMap( ConstDataMap const & src ) : data_map_( src.data_map_ ) {}

ConstDataMap::~ConstDataMap() = default;

/// @details merge-sort style iteration over the elements of both data-maps. This
/// implementation may not be necessary as the STL library may already implement
/// map's assignment operator this way.
ConstDataMap &
ConstDataMap::operator = ( ConstDataMap const & rhs )
{
	if ( this != &rhs ) {
		auto this_cat_iter = data_map_.begin();
		auto rhs_cat_iter = rhs.data_map_.begin();
		while ( this_cat_iter != data_map_.end() && rhs_cat_iter != rhs.data_map_.end() ) {
			if ( this_cat_iter->first == rhs_cat_iter->first ) {
				NamedConstObjectMap & this_ncom( this_cat_iter->second );
				NamedConstObjectMap const & rhs_ncom( rhs_cat_iter->second );

				auto this_ncom_iter = this_ncom.begin();
				auto rhs_ncom_iter = rhs_ncom.begin();
				while ( this_ncom_iter != this_ncom.end() && rhs_ncom_iter != rhs_ncom.end() ) {
					if ( this_ncom_iter->first == rhs_ncom_iter->first ) {
						if ( this_ncom_iter->second.get() != rhs_ncom_iter->second.get() ) {
							// pointer comparison -- no need to copy the pointer if they are already identical
							this_ncom_iter->second = rhs_ncom_iter->second;
						}
						++this_ncom_iter;
						++rhs_ncom_iter;
					} else if ( this_ncom_iter->first < rhs_ncom_iter->first ) {
						auto next_this_ncom_iter( this_ncom_iter );
						++next_this_ncom_iter;
						this_ncom.erase( this_ncom_iter );
						this_ncom_iter = next_this_ncom_iter;
					} else {
						this_ncom[ rhs_ncom_iter->first ] = rhs_ncom_iter->second;
						++rhs_ncom_iter;
					}
				}

				// delete all extra contents from this_ncom
				while ( this_ncom_iter != this_ncom.end() ) {
					auto next_this_ncom_iter( this_ncom_iter );
					++next_this_ncom_iter;
					this_ncom.erase( this_ncom_iter );
					this_ncom_iter = next_this_ncom_iter;
				}

				// copy all extra contents from rhs_ncom
				while ( rhs_ncom_iter != rhs_ncom.end() ) {
					this_ncom[ rhs_ncom_iter->first ] = rhs_ncom_iter->second;
					++rhs_ncom_iter;
				}

				++this_cat_iter;
				++rhs_cat_iter;
			} else if ( this_cat_iter->first < rhs_cat_iter->first ) {
				auto next_this_cat_iter( this_cat_iter );
				++next_this_cat_iter;
				data_map_.erase( this_cat_iter );
				this_cat_iter = next_this_cat_iter;
			} else {
				data_map_[ rhs_cat_iter->first ] = rhs_cat_iter->second;
				++rhs_cat_iter;
			}
		}

		// delete extra elements from data_map_
		while ( this_cat_iter != data_map_.end() ) {
			auto next_this_cat_iter( this_cat_iter );
			++next_this_cat_iter;
			data_map_.erase( this_cat_iter );
			this_cat_iter = next_this_cat_iter;
		}

		// copy extra elements from rhs.data_map_
		while ( rhs_cat_iter != rhs.data_map_.end() ) {
			data_map_[ rhs_cat_iter->first ] = rhs_cat_iter->second;
			++rhs_cat_iter;
		}
	}
	return *this;
}

/// @details merge-sort style iteration over the elements of both data-maps.
bool
ConstDataMap::operator == ( ConstDataMap const & rhs ) const
{
	return data_map_ == rhs.data_map_;

	/*const_iterator this_cat_iter = data_map_.begin();
	const_iterator rhs_cat_iter = rhs.data_map_.begin();
	while ( this_cat_iter != data_map_.end() && rhs_cat_iter != rhs.data_map_.end() ) {
	if ( this_cat_iter->first != rhs_cat_iter->first ) return false;

	NamedConstObjectMap const & this_ncom( this_cat_iter->second );
	NamedConstObjectMap const & rhs_ncom( rhs_cat_iter->second );

	NamedConstObjectMap::const_iterator this_ncom_iter = this_ncom.begin();
	NamedConstObjectMap::const_iterator rhs_ncom_iter = rhs_ncom.begin();
	while ( this_ncom_iter != this_ncom.end() && rhs_ncom_iter != rhs_ncom.end() ) {
	if ( this_ncom_iter->first != rhs_ncom_iter->first ) return false;
	if ( this_ncom_iter->second.get() != rhs_ncom_iter->second.get() ) return false;
	++this_ncom_iter;
	++rhs_ncom_iter;
	}
	if ( this_ncom_iter != this_ncom.end() || rhs_ncom_iter != rhs_ncom.end() ) return false;

	++this_cat_iter;
	++rhs_cat_iter;
	}

	return this_cat_iter == data_map_.end() && rhs_cat_iter == rhs.data_map_.end();*/
}

ConstDataMap::iterator
ConstDataMap::begin() { return data_map_.begin(); }

ConstDataMap::iterator
ConstDataMap::end() { return data_map_.end(); }

ConstDataMap::const_iterator
ConstDataMap::begin() const { return data_map_.begin(); }

ConstDataMap::const_iterator
ConstDataMap::end() const { return data_map_.end(); }

void
ConstDataMap::add( std::string const & type, std::string const & name, utility::pointer::ReferenceCountCOP const op ){
	data_map_[ type ][ name ] = op;
}

bool
ConstDataMap::has( std::string const & type ) const {
	return data_map_.find( type ) != data_map_.end();
}

bool
ConstDataMap::has( std::string const & type, std::string const & name ) const {
	auto it = data_map_.find( type );
	if ( it == data_map_.end() ) return false;

	auto it2 = it->second.find( name );
	if ( it2 == it->second.end() ) return false;

	return true;
}

ConstDataMap::NamedConstObjectMap &
ConstDataMap::operator []( std::string const & type ) {
	return data_map_[ type ];
}

platform::Size
ConstDataMap::size() const {
	platform::Size count = 0;
	 for ( auto const & iter : data_map_ ) {
		count += iter.second.size();
	}
	return count;
}

} // datacache
} // basic

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
basic::datacache::ConstDataMap::save( Archive & arc ) const {
	arc( CEREAL_NVP( data_map_ ) ); // CategorizedConstObjectMap
}

/// @brief Deserialization method that first deserializes into a
/// map of non-constant owning pointers and then copies elements out
/// of this map into the constant-owning-pointer member variable.
template< class Archive >
void
basic::datacache::ConstDataMap::load( Archive & arc ) {
	typedef std::map< std::string, utility::pointer::ReferenceCountOP > NamedObjectMap;
	typedef std::map< std::string, NamedObjectMap > CategorizedObjectMap;

	CategorizedObjectMap local_data_map;
	arc( local_data_map );
	for ( CategorizedObjectMap::const_iterator iter = local_data_map.begin();
			iter != local_data_map.end(); ++iter ) {
		for ( NamedObjectMap::const_iterator inner_iter = iter->second.begin();
				inner_iter != iter->second.end(); ++inner_iter ) {
			data_map_[ iter->first ][ inner_iter->first ] = inner_iter->second;
		}
	}
}

SAVE_AND_LOAD_SERIALIZABLE( basic::datacache::ConstDataMap );
CEREAL_REGISTER_TYPE( basic::datacache::ConstDataMap )

CEREAL_REGISTER_DYNAMIC_INIT( basic_datacache_ConstDataMap )
#endif // SERIALIZATION
