// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/MinimizationData.cc
/// @brief  A container class for use by certain EnergyMethods during derivative and
//          score function evaluation during minimization routines.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit headers
#include <core/scoring/MinimizationData.hh>

// Project headers
#include <basic/datacache/CacheableData.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

ResSingleMinimizationData::ResSingleMinimizationData() : data_cache_( n_min_single_data ) {}
ResSingleMinimizationData::~ResSingleMinimizationData() {}
ResSingleMinimizationData::ResSingleMinimizationData( ResSingleMinimizationData const & other ) :
	utility::pointer::ReferenceCount(),
	data_cache_( n_min_single_data )
{
	for ( Size ii = 1; ii <= n_min_single_data; ++ii ) {
		data_cache_[ ii ] = other.data_cache_[ ii ] ? other.data_cache_[ ii ]->clone() : CacheableDataOP( 0 );
	}
}

ResSingleMinimizationData &
ResSingleMinimizationData::operator = ( ResSingleMinimizationData const & rhs ) {
	if ( this != & rhs ) {
		for ( Size ii = 1; ii <= n_min_single_data; ++ii ) {
			data_cache_[ ii ] = rhs.data_cache_[ ii ] ? rhs.data_cache_[ ii ]->clone() : CacheableDataOP( 0 );
		}
	}
	return *this;
}

void ResSingleMinimizationData::set_data( min_single_data index, CacheableDataOP data ) {
	data_cache_[ index ] = data;
}

ResSingleMinimizationData::CacheableDataOP
ResSingleMinimizationData::get_data( min_single_data index )
{
	return data_cache_[ index ];
}

ResSingleMinimizationData::CacheableDataCOP
ResSingleMinimizationData::get_data( min_single_data index ) const
{
	return data_cache_[ index ];
}


ResPairMinimizationData::ResPairMinimizationData() : data_cache_( n_min_pair_data ) {}
ResPairMinimizationData::~ResPairMinimizationData() {}
ResPairMinimizationData::ResPairMinimizationData( ResPairMinimizationData const & other ) :
	utility::pointer::ReferenceCount(),
	data_cache_( n_min_pair_data )
{
	for ( Size ii = 1; ii <= n_min_pair_data; ++ii ) {
		data_cache_[ ii ] = other.data_cache_[ ii ] ? other.data_cache_[ ii ]->clone() : CacheableDataOP( 0 );
	}
}

ResPairMinimizationData &
ResPairMinimizationData::operator = ( ResPairMinimizationData const & rhs )
{
	if ( this != & rhs ) {
		for ( Size ii = 1; ii <= n_min_pair_data; ++ii ) {
			data_cache_[ ii ] = rhs.data_cache_[ ii ] ? rhs.data_cache_[ ii ]->clone() : CacheableDataOP( 0 );
		}
	}
	return *this;
}

void ResPairMinimizationData::set_data( min_pair_data index, CacheableDataOP data )
{
	data_cache_[ index ] = data;
}

ResPairMinimizationData::CacheableDataOP
ResPairMinimizationData::get_data( min_pair_data index )
{
	return data_cache_[ index ];
}

ResPairMinimizationData::CacheableDataCOP
ResPairMinimizationData::get_data( min_pair_data index ) const
{
	return data_cache_[ index ];
}


}
}

