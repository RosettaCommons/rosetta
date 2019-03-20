// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/geometry/hashing/MinimalClashHash.cc
/// @brief  Minimal hashing class to find clashes as fast as possible
/// @details A rather limiting design choice of this class is to assume
///          that all atoms have the same radius. xyzStripeHash is
///          the accurate version of this class.
/// @author Brian Coventry (bcov@uw.edu)

#include <numeric/geometry/hashing/MinimalClashHash.hh>

#include <numeric/VoxelArray.hh>

#include <utility/pointer/memory.hh>

#include <ObjexxFCL/format.hh>


namespace numeric {
namespace geometry {
namespace hashing {

MinimalClashHash::MinimalClashHash(
	float grid_resolution /*= 0.0*/,
	float ball_radius /*= 2.0f*/,
	utility::vector1<Ball> const & balls/*=utility::vector1<Ball>()*/
) :
	grid_resolution_( grid_resolution ),
	ball_radius_( ball_radius )
{
	init( balls );
}

void
MinimalClashHash::init( utility::vector1<Ball> const & balls ) {

	// Find the edges of the box

	Vec ub{ -1e9, -1e9, -1e9};
	Vec lb{  1e9,  1e9,  1e9};

	for ( Ball const & ball : balls ) {
		Vec const & xyz = ball.xyz();

		for ( int i = 0; i < 3; i++ ) {
			ub[ i ] = std::max<float>( ub[ i ], xyz[ i ] );
			lb[ i ] = std::min<float>( lb[ i ], xyz[ i ] );
		}
	}

	// Pad the edges of the box
	for ( int i = 0; i < 3; i++ ) {
		ub[ i ] += ball_radius_ * 2 + grid_resolution_ * 2; // should theoretically be ok without doing
		lb[ i ] -= ball_radius_ * 2 + grid_resolution_ * 2; // double resolution. But I don't trust it.
	}

	// Make the voxel array at the right resolution
	Vec cs{ grid_resolution_, grid_resolution_, grid_resolution_ };

	voxel_array_ = utility::pointer::make_shared< VoxelArray >( lb, ub, cs );


	// Fill the voxels around each ball whose center lies within 2 x ball_radius

	// This is how many steps we must take on either side to make sure we cover all voxels
	int double_radius_steps = (int) std::ceil( ball_radius_ * 2 / grid_resolution_ ) + 1;

	float double_radius_sq = ball_radius_ * ball_radius_ * 4;

	typedef typename VoxelArray::Bounds Bounds;

	for ( Ball const & ball : balls ) {

		Vec const & xyz = ball.xyz();

		Bounds centered_xyz = voxel_array_->indices_to_center( voxel_array_->floats_to_index( xyz ) );

		Vec worker { 0, 0, 0};
		for ( int ix = -double_radius_steps; ix <= double_radius_steps; ix++ ) {
			worker[0] = ix * grid_resolution_ + centered_xyz[0];

			for ( int iy = -double_radius_steps; iy <= double_radius_steps; iy++ ) {
				worker[1] = iy * grid_resolution_ + centered_xyz[1];

				for ( int iz = -double_radius_steps; iz <= double_radius_steps; iz++ ) {
					worker[2] = iz * grid_resolution_ + centered_xyz[2];

					if ( xyz.distance_squared( worker ) > double_radius_sq ) continue;

					(*voxel_array_)[ worker ] = true;
				}
			}
		}
	}
}


bool
MinimalClashHash::clash( Vec const & v ) const {
	return voxel_array_->at( v );
}
bool
MinimalClashHash::clash_ball( Ball const & b ) const {
	return clash( b.xyz() );
}
Size
MinimalClashHash::clash_check_balls( utility::vector1<Ball> const & bs, Size clash_limit /*= 100000000 */) const {
	Size clashes = 0;
	for ( Ball const & ball : bs ) {
		if ( clash_ball( ball ) ) clashes++;
		if ( clashes >= clash_limit ) break;
	}
	return clashes;
}


} // namespace hashing
} // namespace geometry
} // namespace numeric
