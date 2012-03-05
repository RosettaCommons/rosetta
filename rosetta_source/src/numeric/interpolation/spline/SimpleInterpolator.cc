// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/interpolation/spline/SimpleInterpolator.cc
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler
///

#include <numeric/interpolation/spline/SimpleInterpolator.hh>
#include <numeric/interpolation/spline/spline_functions.hh>

namespace numeric {
namespace interpolation {
namespace spline {

SimpleInterpolator::SimpleInterpolator(
	 utility::vector1<Real> const & x,
	 utility::vector1<Real> const & y,
	 Real lbdy,
	 Real ubdy
) :
	Interpolator(),
	x_(x),
	y_(y),
	ddy_()
{
	ddy_ = spline_second_derivative(x_,y_,lbdy,ubdy);
}

void
SimpleInterpolator::interpolate( Real x, Real & y, Real & dy ) {
	if(has_lb_function() && x < get_lb_function_cutoff())
	{
		return compute_lb_function_solution(x,y);
	}
	if(has_ub_function() && x > get_ub_function_cutoff())
	{
		return compute_ub_function_solution(x,y);
	}
	return spline_interpolate(x_,y_,ddy_,x,y,dy);
}

} // end namespace spline
} // end namespace interpolation
} // end namespace numeric
