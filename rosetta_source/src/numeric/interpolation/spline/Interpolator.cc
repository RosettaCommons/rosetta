// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/interpolation/spline/Interpolator.cc
/// @author Sam DeLuca


#include <numeric/interpolation/spline/Interpolator.hh>

namespace numeric {
namespace interpolation {
namespace spline {

Interpolator::Interpolator() :
	has_lb_function_(false),
	has_ub_function_(false),
	lb_cutoff_(0.0),
	ub_cutoff_(0.0),
	lb_slope_(0.0),
	ub_slope_(0.0),
	lb_intercept_(0.0),
	ub_intercept_(0.0)

{

}

void Interpolator::set_lb_function(Real const &  lb, Real const & slope, Real const & intercept)
{
	lb_cutoff_ = lb;
	lb_slope_ = slope;
	lb_intercept_ = intercept;
	has_lb_function_ = true;
}

void Interpolator::set_ub_function(Real const & ub, Real const & slope,Real const & intercept)
{
	ub_cutoff_ = ub;
	ub_slope_ = slope;
	ub_intercept_ = intercept;
	has_ub_function_ = true;
}

bool Interpolator::has_lb_function() const
{
	return has_lb_function_;
}

bool Interpolator::has_ub_function() const
{
	return has_ub_function_;
}

Real Interpolator::get_lb_function_cutoff() const
{
	return lb_cutoff_;
}

Real Interpolator::get_ub_function_cutoff() const
{
	return ub_cutoff_;
}

void Interpolator::compute_lb_function_solution(Real x, Real & y) const
{
	y = lb_slope_*x+lb_intercept_;
}

void Interpolator::compute_ub_function_solution(Real x, Real & y) const
{
	y = ub_slope_*x+ub_intercept_;
}

}
}
}
