// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableData.cc
/// @brief
/// @author Phil Bradley

// Unit headers
#include <basic/datacache/CacheableData.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#endif

namespace basic {
namespace datacache {

CacheableData::~CacheableData() {}

}
}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void basic::datacache::CacheableData::save( Archive & ) const {}

/// @brief Automatically generated deserialization method
template< class Archive >
void basic::datacache::CacheableData::load( Archive & ) {}

SAVE_AND_LOAD_SERIALIZABLE( basic::datacache::CacheableData );
CEREAL_REGISTER_DYNAMIC_INIT( basic_datacache_CacheableData )

#endif // SERIALIZATION
