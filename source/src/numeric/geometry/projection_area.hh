// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
// //
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
// /// @file    numeric/geometry/projection_area.hh
// /// @brief   Function to get an approximate projection area.
// /// @author  SM Baargeen Alam Turzo
// /// @author  Turzo <turzo.1@osu.edu>
//

#ifndef INCLUDED_numeric_geometry_projection_area_HH
#define INCLUDED_numeric_geometry_projection_area_HH

// Numeric headers
#include <numeric/types.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/MathMatrix.fwd.hh>

// Utility header
#include <utility/vector1.fwd.hh>

/// @brief Return the 2d area of a projection, need x and y (or z) coordinates as well as their respective
/// element radius and probe radius, probe radius can be set to zero if only the projection of protein
/// in vaccum is desired.
namespace numeric {
namespace geometry {
Real projection_area(utility::vector1 < Real > const &xs, utility::vector1 < Real > const &ys, utility::vector1 < Real > const &elements_rad, Real const probe_radius);
}
}
#endif  // INCLUDED_numeric_geometry_projection_area_HH
