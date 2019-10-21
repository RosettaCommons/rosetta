// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableResRotPairFloatMap.hh
/// @brief
/// @author Brian Coventry

// unit headers
#include <basic/datacache/CacheableResRotPairFloatMap.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

#include <numeric/MathMatrix.srlz.hh>

// Cereal headers
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {

CacheableResRotPairFloatMap::CacheableResRotPairFloatMap() : CacheableData() {}

CacheableResRotPairFloatMap::~CacheableResRotPairFloatMap() = default;

CacheableDataOP
CacheableResRotPairFloatMap::clone() const {
	return CacheableDataOP( new CacheableResRotPairFloatMap(*this) );
}

CacheableResRotPairFloatMapOP
CacheableResRotPairFloatMap::shared_from_this() {
	return utility::pointer::static_pointer_cast<CacheableResRotPairFloatMap>( CacheableData::shared_from_this() );
}

std::unordered_map< ResRotPair, float, ResRotPairHasher > &
CacheableResRotPairFloatMap::map(){
	return map_;
}

std::unordered_map< ResRotPair, float, ResRotPairHasher > const &
CacheableResRotPairFloatMap::map() const {
	return map_;
}


} // namespace datacache
} // namespace basic


#ifdef    SERIALIZATION
/// @brief Automatically generated serialization method
template< class Archive >
void
basic::datacache::ResRotPair::save( Archive & arc ) const {
	arc( CEREAL_NVP( first_res ) );
	arc( CEREAL_NVP( first_rot ) );
	arc( CEREAL_NVP( second_res ) );
	arc( CEREAL_NVP( second_rot ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
basic::datacache::ResRotPair::load( Archive & arc ) {
	arc( first_res );
	arc( first_rot );
	arc( second_res );
	arc( second_rot );
}

/// @brief Automatically generated serialization method
template< class Archive >
void
basic::datacache::CacheableResRotPairFloatMap::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( map_ ) ); // std::unordered_map<uint64_t, MathMatrix<float>>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
basic::datacache::CacheableResRotPairFloatMap::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( map_ ); // std::unordered_map<uint64_t, MathMatrix<float>>
}

SAVE_AND_LOAD_SERIALIZABLE( basic::datacache::CacheableResRotPairFloatMap );
CEREAL_REGISTER_TYPE( basic::datacache::CacheableResRotPairFloatMap )

CEREAL_REGISTER_DYNAMIC_INIT( basic_datacache_CacheableResRotPairFloatMap )
#endif // SERIALIZATION
