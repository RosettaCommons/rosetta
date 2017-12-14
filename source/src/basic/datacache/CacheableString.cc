// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/datacache/CacheableString.hh
/// @brief
/// @author Phil Bradley


// unit headers
#include <basic/datacache/CacheableString.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace basic {
namespace datacache {

CacheableString::CacheableString( std::string const & str ) : CacheableData(), str_(str) {}
CacheableString::~CacheableString() = default;
CacheableDataOP CacheableString::clone() const { return CacheableDataOP( new CacheableString(*this) ); }
std::string const & CacheableString::str() const { return str_; }

CacheableStringOP CacheableString::shared_from_this() { return utility::pointer::static_pointer_cast<CacheableString>( CacheableData::shared_from_this() ); }

} // namespace datacache
} // namespace basic


#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
basic::datacache::CacheableString::CacheableString() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
basic::datacache::CacheableString::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( str_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
basic::datacache::CacheableString::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( str_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( basic::datacache::CacheableString );
CEREAL_REGISTER_TYPE( basic::datacache::CacheableString )

CEREAL_REGISTER_DYNAMIC_INIT( basic_datacache_CacheableString )
#endif // SERIALIZATION
