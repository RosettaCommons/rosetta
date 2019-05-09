// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/linear_algebra/minimum_bounding_ellipse.hh
/// @brief Given a set of points, calculate the minimum bounding ellipse
/// @author Rebecca Alford (ralford3@jhu.edu)

#ifndef INCLUDED_numeric_linear_algebra_minimum_bounding_ellipse_HH
#define INCLUDED_numeric_linear_algebra_minimum_bounding_ellipse_HH

// Project Headers
#include <numeric/MathMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/types.hh>

#include <numeric/linear_algebra/EllipseParameters.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <cstdlib>

namespace numeric {
namespace linear_algebra {


/// @brief Use the Khachiyan Algorithm to compute the minimum volume enclosing ellipsoid given a set of (x,y) data points
EllipseParametersOP
minimum_bounding_ellipse(
	utility::vector1< xyzVector< Real > > points,
	Real tolerance,
	Size max_iterations=50
);

/// @brief Calculate the sum-of-square differences between
/// values stored in two vector1 objects
Real
sum_of_square_differences(
	MathMatrix< Real > old_u,
	MathMatrix< Real > new_u
);

/// @brief Calculate the transpose of a non-square MathMatrix
/// and return the result as a new MathMatrix
MathMatrix< Real >
non_square_transpose(
	MathMatrix< Real > matrix_in
);

/// @brief Check whether a given test point lies within an ellipse
bool
point_in_ellipse(
	xyzVector< Real > p,
	Real const h,
	Real const k,
	Real const a,
	Real const b,
	MathMatrix< Real > rotation
);

} // ns linear_algebra
} // ns numeric

#endif // INCLUDED_mumeric_linear_algebra_minimum_bounding_ellipse_HH
