// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableStringIntegerMap.hh
/// @brief
/// @author Phil Bradley, Jared Adolf-Bryfogle (StringInt)

// unit headers
#include <basic/datacache/CacheableStringIntegerMap.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {

CacheableStringIntegerMap::CacheableStringIntegerMap() : CacheableData() {}

CacheableStringIntegerMap::~CacheableStringIntegerMap() = default;

CacheableDataOP
CacheableStringIntegerMap::clone() const { return CacheableDataOP( new CacheableStringIntegerMap(*this) ); }

CacheableStringIntegerMapOP
CacheableStringIntegerMap::shared_from_this() { return utility::pointer::static_pointer_cast<CacheableStringIntegerMap>( CacheableData::shared_from_this() ); }

std::map< std::string, int > &
CacheableStringIntegerMap::map(){ return map_; }

std::map< std::string, int > const &
CacheableStringIntegerMap::map() const { return map_; }


} // namespace datacache
} // namespace basic


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
basic::datacache::CacheableStringIntegerMap::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( map_ ) ); // std::map<std::string, int>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
basic::datacache::CacheableStringIntegerMap::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( map_ ); // std::map<std::string, int>
}

SAVE_AND_LOAD_SERIALIZABLE( basic::datacache::CacheableStringIntegerMap );
CEREAL_REGISTER_TYPE( basic::datacache::CacheableStringIntegerMap )

CEREAL_REGISTER_DYNAMIC_INIT( basic_datacache_CacheableStringIntegerMap )
#endif // SERIALIZATION
