// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/SixDHasher.cc
/// @brief  Implementation of 6D hasher classes
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini


// Unit headers
#include <protocols/match/SixDHasher.hh>


namespace protocols {
namespace match {


SixDCoordinateBinner::SixDCoordinateBinner()
{}

SixDCoordinateBinner::SixDCoordinateBinner(
	BoundingBox const & bounding_box,
	Size3 const & euler_offsets,
	utility::fixedsizearray1< Real, 6 > bin_widths
) :
	bounding_box_( bounding_box ),
	bin_widths_( bin_widths )
{
	Vector span = bounding_box_.upper() - bounding_box_.lower();
	Vector new_upper = bounding_box_.upper();

	for ( Size ii = 1; ii <= 3; ++ii ) {
		dimsizes_[ ii ] = static_cast< Size > ( span( ii ) / bin_widths_[ ii ] );
		if ( dimsizes_[ ii ] * bin_widths_[ ii ] < span( ii ) ) {
			dimsizes_[ ii ] += 1;
			new_upper( ii ) = bounding_box_.lower()( ii ) + dimsizes_[ ii ] * bin_widths_[ ii ];
		}
	}
	bounding_box_.set_upper( new_upper );

	for ( Size ii = 4; ii <= 5; ++ii ) {
		dimsizes_[ ii ] = static_cast< Size > ( 360.0 / bin_widths_[ ii ] );
	}
	dimsizes_[ 6 ] = static_cast< Size > ( 180.0 / bin_widths_[ 6 ] );
	if( dimsizes_[6] == 0 ) dimsizes_[6] = 1;

	/// Add an extra bin so that values near 180 ( wi binwidth/2) and near 0 ( wi binwidth/2 ) can be joined.
	if ( euler_offsets[ 3 ] != 0 ) {
		++dimsizes_[ 6 ];
	}

	dimprods_[ 6 ] = 1;
	for ( Size ii = 5; ii >= 1; --ii ) {
		dimprods_[ ii ] = dimprods_[ ii + 1 ] * dimsizes_[ ii + 1 ];
	}

	for ( Size ii = 1; ii <= 3; ++ii ) {
		if ( euler_offsets[ ii ] != 0 ) {
			euler_offsets_[ ii ] = bin_widths_[ ii + 3 ] / 2;
		} else {
			euler_offsets_[ ii ] = 0.0;
		}
	}
	for ( Size ii = 1; ii <= 6; ++ii ) { halfbin_widths_[ ii ] = bin_widths_[ ii ] / 2; }
}


/// @details When floating point comparison breaks down, it is possible to have a point
/// hash outside of the bounding volume:  359.999999 / 10.0 -> 36 instead of 35.
/// For this reason, it's important to mod the resulting bin index by the number of bins.
Bin6D
SixDCoordinateBinner::bin6( Real6 const & values ) const {

	core::Vector xyzcoord( values[ 1 ], values[ 2 ], values[ 3 ] );

	assert( bounding_box_.contains( xyzcoord ) );
	assert( values[ 4 ] >= 0.0 && values[ 4 ] < 360.0 );
	assert( values[ 5 ] >= 0.0 && values[ 5 ] < 360.0 );
	assert( values[ 6 ] >= 0 && values[ 6 ] <= 180.0 );

	core::Vector from_corner = xyzcoord - bounding_box_.lower();
	Bin6D bins;

	bins[ 1 ] = static_cast< Size > ( from_corner.x() / bin_widths_[ 1 ] );
	if ( bins[ 1 ] == dimsizes_[ 1 ] ) --bins[ 1 ];

	bins[ 2 ] = static_cast< Size > ( from_corner.y() / bin_widths_[ 2 ] );
	if ( bins[ 2 ] == dimsizes_[ 2 ] ) --bins[ 2 ];

	bins[ 3 ] = static_cast< Size > ( from_corner.z() / bin_widths_[ 3 ] );
	if ( bins[ 3 ] == dimsizes_[ 3 ] ) --bins[ 3 ];

	Real3 euler = wrap_euler_angles( values );

	bins[ 4 ] = static_cast< Size > ( euler[ 1 ] / bin_widths_[ 4 ] ) % dimsizes_[ 4 ];
	bins[ 5 ] = static_cast< Size > ( euler[ 2 ] / bin_widths_[ 5 ] ) % dimsizes_[ 5 ];
	bins[ 6 ] = static_cast< Size > ( euler[ 3 ] / bin_widths_[ 6 ] ) % dimsizes_[ 6 ];

	return bins;
}


/// @details
Bin6D
SixDCoordinateBinner::halfbin6( Real6 const & values ) const
{
	Bin6D bin = bin6( values );

	core::Vector xyzcoord( values[ 1 ], values[ 2 ], values[ 3 ] );
	core::Vector from_corner = xyzcoord - bounding_box_.lower();

	core::Vector r3lower_corner = bounding_box_.lower() +
		Vector( bin[1]*bin_widths_[1], bin[2]*bin_widths_[2], bin[3]*bin_widths_[3] );
	Vector remnant = xyzcoord - r3lower_corner;

	Bin6D halfbin( 0 );
	for ( Size ii = 1; ii <= 3; ++ii ) {
		if ( remnant[ ii-1 ] >= halfbin_widths_[ ii ]) {
			halfbin[ii] = 1;
		}
	}

	Real3 euler = wrap_euler_angles( values );
	for ( Size ii = 1; ii <= 3; ++ii ) {
		Real lower_boundary = bin[ii+3]*bin_widths_[ii+3];
		Real dist = euler[ii] - lower_boundary;
		//std::cout << "Euler " << ii << " " << values[ii+3] <<  " " << euler[ ii ] << " " << lower_boundary << " " << dist << " " << euler_offsets_[ii] << std::endl;
		if ( dist >= halfbin_widths_[ii+3] ) {
			halfbin[ii+3] = 1;
		}
	}
	return halfbin;
}

Real6
SixDCoordinateBinner::bin_center_point( Bin6D const & bin ) const
{
	Real6 center;
	for ( Size ii = 1; ii <= 3; ++ii ) center[ ii ] = bounding_box_.lower()( ii ) + bin[ ii ] * bin_widths_[ ii ] + halfbin_widths_[ ii ];
	for ( Size ii = 1; ii <= 3; ++ii ) center[ ii + 3 ] = euler_offsets_[ ii ] + bin[ ii + 3 ] * bin_widths_[ ii + 3 ] + halfbin_widths_[ ii + 3 ];
	return center;
}


Real3
SixDCoordinateBinner::wrap_euler_angles( Real6 const & values ) const
{
	Real3 euler;

	if ( values[ 6 ] > euler_offsets_[ 3 ] ) {
		euler[ 3 ] = values[ 6 ] + euler_offsets_[ 3 ];

	} else {
		euler[ 3 ] = values[ 6 ];
	}

	if (( values[ 6 ] < euler_offsets_[ 3 ] || values[ 6 ] >= 180.0 - euler_offsets_[ 3 ] ) &&
			( values[ 4 ] - euler_offsets_[ 1 ] > 180.0 )) {
		/// WRAP if phi > 180 to the region of negative theta rotation.
		/// The idea is that, when we're wrapping theta, half of the points have to stay fixed so they end up in the same
		/// bin: if all points wrapped, then none would land in the same bin.
		euler[ 1 ] = values[ 4 ] - euler_offsets_[ 1 ] - 180.0;
		euler[ 2 ] = values[ 5 ] - euler_offsets_[ 2 ] - 180.0;
		if ( euler[ 2 ] < 0 ) euler[ 2 ] += 360;
	} else {
		/// Leave phi/psi in their usual positions
		euler[ 1 ] = values[ 4 ] - euler_offsets_[ 1 ];
		if ( euler[ 1 ] < 0 ) euler[ 1 ] += 360;

		euler[ 2 ] = values[ 5 ] - euler_offsets_[ 2 ];
		if ( euler[ 2 ] < 0 ) euler[ 2 ] += 360;
	}
	return euler;
}




}
}

