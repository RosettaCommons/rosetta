// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/interpolation/spline/SimpleInterpolator.cc
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler


#include <numeric/interpolation/spline/SimpleInterpolator.hh>
#include <numeric/interpolation/spline/spline_functions.hh>

#include <utility/tools/make_vector.hh>

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

SimpleInterpolator::SimpleInterpolator() :Interpolator()
{

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

utility::json_spirit::Value SimpleInterpolator::serialize()
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	std::vector<Value> x_values,y_values,ddy_values;

	for(utility::vector1<Real>::iterator it = x_.begin(); it != x_.end();++it)
	{
		x_values.push_back(Value(*it));
	}

	for(utility::vector1<Real>::iterator it = y_.begin(); it != y_.end();++it)
	{
		y_values.push_back(Value(*it));
	}

	for(utility::vector1<Real>::iterator it = ddy_.begin(); it != ddy_.end();++it)
	{
		ddy_values.push_back(Value(*it));
	}

	Pair x_data("xdata",x_values);
	Pair y_data("ydata",y_values);
	Pair ddy_data("ddydata",ddy_values);

	Pair base_data("base_data",Interpolator::serialize());

	return Value(utility::tools::make_vector(x_data,y_data,ddy_data,base_data));

}

void SimpleInterpolator::deserialize(utility::json_spirit::mObject data)
{
	utility::json_spirit::mArray x_data(data["xdata"].get_array());
	utility::json_spirit::mArray y_data(data["ydata"].get_array());
	utility::json_spirit::mArray ddy_data(data["ddydata"].get_array());

	x_.clear();
	y_.clear();
	ddy_.clear();

	for(utility::json_spirit::mArray::iterator it = x_data.begin();it != x_data.end();++it)
	{
		x_.push_back(it->get_real());
	}

	for(utility::json_spirit::mArray::iterator it = y_data.begin();it != y_data.end();++it)
	{
		y_.push_back(it->get_real());
	}

	for(utility::json_spirit::mArray::iterator it = ddy_data.begin();it != ddy_data.end();++it)
	{
		ddy_.push_back(it->get_real());
	}

	Interpolator::deserialize(data["base_data"].get_obj());

}


} // end namespace spline
} // end namespace interpolation
} // end namespace numeric
