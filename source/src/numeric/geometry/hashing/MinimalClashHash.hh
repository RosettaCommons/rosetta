// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/geometry/hashing/MinimalClashHash.hh
/// @brief  Minimal hashing class to find clashes as fast as possible
/// @details A rather limiting design choice of this class is to assume
///          that all atoms have the same radius. xyzStripeHash is
///          the accurate version of this class.
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_numeric_geometry_hashing_MinimalClashHash_hh
#define INCLUDED_numeric_geometry_hashing_MinimalClashHash_hh

#include <numeric/geometry/hashing/xyzStripeHash.hh> // need Ball
#include <numeric/VoxelArray.fwd.hh>  // Try to avoid including the .hh here
#include <numeric/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

namespace numeric {
namespace geometry {
namespace hashing {

class MinimalClashHash : public utility::pointer::ReferenceCount {

	typedef numeric::xyzVector<float> Vec;
	typedef numeric::VoxelArray<float, bool> VoxelArray;
	typedef utility::pointer::shared_ptr< VoxelArray > VoxelArrayOP;

public:
	MinimalClashHash(
		float grid_resolution = 0.0f,
		float ball_radius = 2.0f,
		utility::vector1<Ball> const & balls=utility::vector1<Ball>()
	);

	void init( utility::vector1<Ball> const & balls );


	bool clash( Vec const & v ) const;
	bool clash_ball( Ball const & b ) const;
	Size clash_check_balls( utility::vector1<Ball> const & bs, Size clash_limit = 100000000 ) const;


private:
	float grid_resolution_;
	float ball_radius_;

	VoxelArrayOP voxel_array_;
};


} // namespace hashing
} // namespace geometry
} // namespace numeric

#endif
