// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/pack/rotamer_set/WaterPackingInfo.cc
/// @brief
/// @author


#include <core/types.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.hh>

// utility headers

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace rotamer_set {

WaterPackingInfo::WaterPackingInfo() = default;

WaterPackingInfo::WaterPackingInfo( WaterPackingInfo const & src ):
	CacheableData(),
	data_( src.data_ )
{
	for ( Size i=1; i<= data_.size(); ++i ) {
		if ( data_[ i ] ) data_[i] = WaterAnchorInfoOP( new WaterAnchorInfo( *data_[i] ) );
	}
}

basic::datacache::CacheableDataOP
WaterPackingInfo::clone() const {
	return basic::datacache::CacheableDataOP( new WaterPackingInfo( *this ) );
}

WaterAnchorInfo &
WaterPackingInfo::operator[] ( Size const seqpos ) {
	if ( seqpos > data_.size() ) data_.resize( seqpos, nullptr );
	if ( data_[seqpos] == nullptr ) {
		data_[seqpos] = WaterAnchorInfoOP( new WaterAnchorInfo() );
	}
	return *( data_[ seqpos ] );
}

WaterAnchorInfo const &
WaterPackingInfo::operator[] ( Size const seqpos ) const {
	debug_assert( seqpos <= data_.size() && data_[ seqpos ] );
	return *( data_[ seqpos ] );
}


void
WaterPackingInfo::clear() {
	data_.clear();
}

} // namespace rotamer_set
} // namespace pack
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::rotamer_set::WaterPackingInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( data_ ) ); // utility::vector1<WaterAnchorInfoOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::rotamer_set::WaterPackingInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( data_ ); // utility::vector1<WaterAnchorInfoOP>
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::rotamer_set::WaterPackingInfo );
CEREAL_REGISTER_TYPE( core::pack::rotamer_set::WaterPackingInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_rotamer_set_WaterPackingInfo )
#endif // SERIALIZATION
