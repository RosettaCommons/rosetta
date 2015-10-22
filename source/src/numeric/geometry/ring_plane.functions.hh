// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    numeric/geometry/ring_plane.functions.hh
/// @brief   Function declarations for ring-plane geometry functions.
/// @author  Brian Weitzner
/// @author  Michael Pacella
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_numeric_geometry_ring_plane_functions_HH
#define INCLUDED_numeric_geometry_ring_plane_functions_HH

// Numeric headers
#include <numeric/types.hh>
#include <numeric/xyzVector.fwd.hh>

// Utility header
#include <utility/vector1.hh>


namespace numeric {
namespace geometry {

/// @brief A zero-length vector to represent a non-plane or a point at the origin.
extern xyzVector< Real > const ZERO_VECTOR;

/// @brief Return the R-squared value, indicating how well one or more points lie in a given plane, if those points were
/// shifted such that their centroid were at the origin.
Real residual_squared_of_points_to_plane(
	utility::vector1< xyzVector< Real > > const & point_coords,
	xyzVector< Real > const & vector_normal_to_plane );

bool are_coplanar( utility::vector1< xyzVector< Real > > const & ring_point_coords );

/// @brief Return the vector normal to the plane of best fit of a ring.
xyzVector< Real > vector_normal_to_ring_plane_of_best_fit(
	utility::vector1< xyzVector< Real > > const & ring_point_coords,
	bool co_planar_check = true );

}  // namespace geometry
}  // namespace numeric

#endif  // INCLUDED_numeric_geometry_ring_plane_functions_HH
