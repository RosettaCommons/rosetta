// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/WriteableCacheableData.hh
/// @brief
/// @author Justin Porter

// unit headers
#include <basic/datacache/WriteableCacheableData.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {

WriteableCacheableData::~WriteableCacheableData() = default;

WriteableCacheableDataOP WriteableCacheableData::shared_from_this() { return utility::pointer::static_pointer_cast<WriteableCacheableData>( CacheableData::shared_from_this() ); }


} // namespace datacache
} // namespace basic


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
basic::datacache::WriteableCacheableData::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
basic::datacache::WriteableCacheableData::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( basic::datacache::WriteableCacheableData );
CEREAL_REGISTER_TYPE( basic::datacache::WriteableCacheableData )

CEREAL_REGISTER_DYNAMIC_INIT( basic_datacache_WriteableCacheableData )
#endif // SERIALIZATION
