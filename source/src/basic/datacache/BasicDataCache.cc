// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/datacache/BasicDataCache.cc
/// @brief  A DataCache storing objects derived from
///         basic::datacache::CacheableData.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <basic/datacache/BasicDataCache.hh>

#include <cstddef>  // for size_t

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal serialization headers
#include <cereal/cereal.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace basic {
namespace datacache {


/// @brief default constructor
BasicDataCache::BasicDataCache() :
	Super()
{}


/// @brief size constructor
/// @param[in] n_types The number of slots for this DataCache.
BasicDataCache::BasicDataCache( size_t const n_slots ) :
	Super( n_slots )
{}


/// @brief copy constructor
BasicDataCache::BasicDataCache( BasicDataCache const & rval ) :
	Super( rval )
{}


/// @brief default destructor
BasicDataCache::~BasicDataCache() {}


/// @brief copy assignment
BasicDataCache & BasicDataCache::operator =( BasicDataCache const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );
	}
	return *this;
}


} // namespace datacache
} // namespace basic

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
basic::datacache::BasicDataCache::save( Archive & arc ) const {
	arc( data() );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
basic::datacache::BasicDataCache::load( Archive & arc ) {
	arc( data() );
}

SAVE_AND_LOAD_SERIALIZABLE( basic::datacache::BasicDataCache );
CEREAL_REGISTER_TYPE( basic::datacache::BasicDataCache )

CEREAL_REGISTER_DYNAMIC_INIT( basic_datacache_BasicDataCache )
#endif // SERIALIZATION
