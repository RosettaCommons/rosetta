// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableUint64MathMatrixFloatMap.hh
/// @brief
/// @author Brian Coventry

// unit headers
#include <basic/datacache/CacheableUint64MathMatrixFloatMap.hh>


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

CacheableUint64MathMatrixFloatMap::CacheableUint64MathMatrixFloatMap() : CacheableData() {}

CacheableUint64MathMatrixFloatMap::~CacheableUint64MathMatrixFloatMap() = default;

CacheableDataOP
CacheableUint64MathMatrixFloatMap::clone() const {
	return CacheableDataOP( new CacheableUint64MathMatrixFloatMap(*this) );
}

CacheableUint64MathMatrixFloatMapOP
CacheableUint64MathMatrixFloatMap::shared_from_this() {
	return utility::pointer::static_pointer_cast<CacheableUint64MathMatrixFloatMap>( CacheableData::shared_from_this() );
}

std::unordered_map< uint64_t, numeric::MathMatrix<float> > &
CacheableUint64MathMatrixFloatMap::map(){
	return map_;
}

std::unordered_map< uint64_t, numeric::MathMatrix<float> > const &
CacheableUint64MathMatrixFloatMap::map() const {
	return map_;
}


} // namespace datacache
} // namespace basic


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
basic::datacache::CacheableUint64MathMatrixFloatMap::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( map_ ) ); // std::unordered_map<uint64_t, MathMatrix<float>>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
basic::datacache::CacheableUint64MathMatrixFloatMap::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( map_ ); // std::unordered_map<uint64_t, MathMatrix<float>>
}

SAVE_AND_LOAD_SERIALIZABLE( basic::datacache::CacheableUint64MathMatrixFloatMap );
CEREAL_REGISTER_TYPE( basic::datacache::CacheableUint64MathMatrixFloatMap )

CEREAL_REGISTER_DYNAMIC_INIT( basic_datacache_CacheableUint64MathMatrixFloatMap )
#endif // SERIALIZATION
