// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    numeric/geometry/ring_plane.functions.hh
/// @brief   Function definitions for ring-plane geometry functions.
/// @author  Brian Weitzner
/// @author  Michael Pacella
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <numeric/geometry/ring_plane.functions.hh>

// Numeric headers
#include <numeric/linear_algebra/singular_value_decomposition.hh>
#include <numeric/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

// Utility header
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>


namespace numeric {
namespace geometry {

// A zero-length vector to represent a non-plane.
xyzVector< Real > const ZERO_VECTOR = xyzVector< Real >( 0.0 );


// Return the R-squared value, indicating how well one or more points lie in a given plane, if those points were
// shifted such that their centroid were at the origin.
/// @details This function takes a list of coordinates of points and measures and sums the residual values of each
/// relative to the plane, which is provided by its normal vector as a parameter to the function.
/// @param   <point_coords>: A list of Cartesian coordinates for one or more points.
/// @param   <vector_normal_to_plane>: The vector describing the plane, which is assumed to pass through the origin.
/// @return  The sum of the squares of the R values of each point relative to the plane. If A) the given points all
/// lie in the same plane and B) that plane is oriented parallel to the given plane, R-squared will be zero.
Real
residual_squared_of_points_to_plane(
	utility::vector1< xyzVector< Real > > const & point_coords,
	xyzVector< Real > const & vector_normal_to_plane )
{
	// First, check that we actually have a plane with which to work.
	if ( vector_normal_to_plane == ZERO_VECTOR ) {
		return 0.0;  // There is no plane to which the points in the ring may be compared!
	}

	// Next, calculate the average point.
	xyzVector< Real > average_point( center_of_mass( point_coords ) );

	// Then, for each point, calculate the displacement vectors from the average point, take its dot product with
	// the normal vector, square it, and add it to the sum.
	Size const n_points( point_coords.size() );
	Real sum_R_squared( 0.0 );
	for ( Size i( 1 ); i <= n_points; ++i ) {
		xyzVector< Real > const disp_vector( point_coords[ i ] - average_point );
		Real const R( vector_normal_to_plane.dot( disp_vector ) );
		Real const R_squared( R * R );
		sum_R_squared += R_squared;
	}

	return sum_R_squared;
}

xyzVector< numeric::Real >
make_plane_from_three_points( xyzVector< Real > const & p1, xyzVector< Real > const & p2, xyzVector< Real > const & p3 )
{
	//define a plane using 3 points
	xyzVector< Real > vec1 = p2 - p1, vec2 = p3 - p1;
	xyzVector< Real > normal_vector = vec1.cross( vec2 );
	return normal_vector.normalize();
}

bool are_coplanar(utility::vector1< xyzVector< Real > > const & ring_point_coords)
{
	// 3 or less points are always coplanar
	if ( ring_point_coords.size() <= 3 ) { return true; }

	// define an origin as the first point
	xyzVector< Real > origin = ring_point_coords[ 1 ];

	// define a plane with the first 3 points
	xyzVector< Real > normal_vector = make_plane_from_three_points( ring_point_coords[ 1 ],
		ring_point_coords[ 2 ],ring_point_coords[ 3 ] );

	// loop through the remaining points, check if vector from p1 to point is orthogonal to the normal vector of the plane
	for ( Size i = 4; i <= ring_point_coords.size(); ++i ) {
		xyzVector< Real > point_to_origin = ring_point_coords[ i ] - origin;
		if ( ( std::abs ( ( point_to_origin.dot( normal_vector ) ) ) > .001) ) { return false; }
	}
	return true;
}

// Return the vector normal to the plane of best fit of a ring.
/// @details This function takes a list of coordinates of the points of a ring of any size, calculates the plane
/// of best fit, and returns the vector normal to that plane.
/// @param   <ring_point_coords>: A list of Cartesian coordinates for the points of a monocyclic ring system in
/// sequence.
/// @return  The normal vector to the plane of best fit or (0.0, 0.0, 0.0) if the list of coordinates given contains less
/// than three points.
xyzVector< Real >
vector_normal_to_ring_plane_of_best_fit( utility::vector1< xyzVector< Real > > const & ring_point_coords, bool co_planar_check )
{
	using utility::vector1;

	// fitting a plane only makes sense if we have 3 or more points
	if ( ring_point_coords.size() < 3 ) { return ZERO_VECTOR; }

	// If we have perfectly coplanar points already (always true for 3 points, check if true for more) our linear system may be
	// ill-conditioned.  In both these cases we can just make a plane from three points in the system instead
	// of least squares regression
	if ( co_planar_check && are_coplanar( ring_point_coords ) ) {
		return make_plane_from_three_points(ring_point_coords[ 1 ], ring_point_coords[ 2 ], ring_point_coords[ 3 ] );
	}

	// subtract out the centroid so that the calculated plane passes through the origin
	xyzVector< Real > const centroid = center_of_mass( ring_point_coords );
	Real const _INITIAL_VALUE_( ( Real ) INT_MAX / 10. );

	Size const m( ring_point_coords.size() ); // length of coordinates vector is same as number of rows
	Size const n( 3 );  // 3 dimensions -- space!

	vector1< vector1< Real > > centered_ring_point_coords( m, vector1< Real >( n,  _INITIAL_VALUE_ ) );
	for ( Size i = 1; i <= ring_point_coords.size(); ++i ) {
		xyzVector< Real > centered_point = ring_point_coords[ i ] - centroid;
		vector1< Real > & centered_point_vector1 = centered_ring_point_coords[ i ];
		centered_point_vector1[ 1 ] = centered_point.x();
		centered_point_vector1[ 2 ] = centered_point.y();
		centered_point_vector1[ 3 ] = centered_point.z();
	}

	vector1< Real > w( n, _INITIAL_VALUE_ );
	vector1< vector1< Real > > v( n, vector1< Real >( n, _INITIAL_VALUE_ ) ); // n x n matrix

	linear_algebra::svdcmp( centered_ring_point_coords, m, n, w, v);

	// index of the smallest singular value
	Size index = arg_min( w );

	//return the singular vector associated with the smallest singular value

	return xyzVector< Real >( v[ 1 ][ index ], v[ 2 ][ index ], v[ 3 ][ index ] );
}


}  // namespace geometry
}  // namespace numeric
