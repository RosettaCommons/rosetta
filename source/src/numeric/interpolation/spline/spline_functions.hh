// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/interpolation/spline_functions.hh
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler


#ifndef INCLUDED_numeric_interpolation_spline_spline_functions_hh
#define INCLUDED_numeric_interpolation_spline_spline_functions_hh

#include <numeric/types.hh>

#include <utility/vector1.hh>


namespace numeric {
namespace interpolation {
namespace spline {


utility::vector1<Real>
spline_second_derivative(
	utility::vector1<Real> const & x,
	utility::vector1<Real> const & y,
	Real yp1,
	Real ypn
);


void
spline_interpolate(
	utility::vector1<Real> const & xa,
	utility::vector1<Real> const & ya,
	utility::vector1<Real> const & y2a,
	Real x, Real & y, Real & dy
);


} // end namespace spline
} // end namespace interpolation
} // end namespace numeric

#endif
