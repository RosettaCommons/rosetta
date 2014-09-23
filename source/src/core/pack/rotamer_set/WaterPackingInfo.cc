// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/pack/rotamer_set/WaterPackingInfo.cc
/// @brief
/// @author


#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/ResidueType.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/pack/rotamer_set/WaterAnchorInfo.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.hh>

// utility headers

namespace core {
namespace pack {
namespace rotamer_set {

WaterPackingInfo::WaterPackingInfo() {}

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
	return new WaterPackingInfo( *this );
}

WaterAnchorInfo &
WaterPackingInfo::operator[] ( Size const seqpos ) {
	if ( seqpos > data_.size() ) data_.resize( seqpos, 0 );
	if ( data_[seqpos] == 0 ) {
		data_[seqpos] = WaterAnchorInfoOP( new WaterAnchorInfo() );
	}
	return *( data_[ seqpos ] );
}

WaterAnchorInfo const &
WaterPackingInfo::operator[] ( Size const seqpos ) const {
	assert( seqpos <= data_.size() && data_[ seqpos ] );
	return *( data_[ seqpos ] );
}


void
WaterPackingInfo::clear() {
	data_.clear();
}

} // namespace rotamer_set
} // namespace pack
} // namespace core
