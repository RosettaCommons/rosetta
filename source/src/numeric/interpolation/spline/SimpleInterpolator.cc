// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/interpolation/spline/SimpleInterpolator.cc
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler


#include <numeric/interpolation/spline/SimpleInterpolator.hh>
#include <numeric/interpolation/spline/spline_functions.hh>

#include <utility/tools/make_vector.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

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
{}

InterpolatorOP
SimpleInterpolator::clone() const {
	return InterpolatorOP( new SimpleInterpolator( *this ) );
}


void
SimpleInterpolator::interpolate( Real x, Real & y, Real & dy ) const {
	if ( has_lb_function() && x < get_lb_function_cutoff() ) {
		return compute_lb_function_solution(x,y);
	}
	if ( has_ub_function() && x > get_ub_function_cutoff() ) {
		return compute_ub_function_solution(x,y);
	}
	return spline_interpolate(x_,y_,ddy_,x,y,dy);
}

utility::json_spirit::Value SimpleInterpolator::serialize() const
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	std::vector<Value> x_values,y_values,ddy_values;

	for ( utility::vector1<Real>::const_iterator it = x_.begin(); it != x_.end(); ++it ) {
		x_values.push_back(Value(*it));
	}

	for ( utility::vector1<Real>::const_iterator it = y_.begin(); it != y_.end(); ++it ) {
		y_values.push_back(Value(*it));
	}

	for ( utility::vector1<Real>::const_iterator it = ddy_.begin(); it != ddy_.end(); ++it ) {
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

	for ( utility::json_spirit::mArray::iterator it = x_data.begin(); it != x_data.end(); ++it ) {
		x_.push_back(it->get_real());
	}

	for ( utility::json_spirit::mArray::iterator it = y_data.begin(); it != y_data.end(); ++it ) {
		y_.push_back(it->get_real());
	}

	for ( utility::json_spirit::mArray::iterator it = ddy_data.begin(); it != ddy_data.end(); ++it ) {
		ddy_.push_back(it->get_real());
	}

	Interpolator::deserialize(data["base_data"].get_obj());

}

bool SimpleInterpolator::operator == ( Interpolator const & other ) const
{
	if ( ! Interpolator::operator==( other ) ) return false;
	SimpleInterpolator const & other_downcast( static_cast< SimpleInterpolator const & > ( other ) );
	if ( x_   != other_downcast.x_   ) return false;
	if ( y_   != other_downcast.y_   ) return false;
	if ( ddy_ != other_downcast.ddy_ ) return false;
	return true;
}

bool SimpleInterpolator::same_type_as_me( Interpolator const & other ) const
{
	return dynamic_cast< SimpleInterpolator const * > (&other);
}


} // end namespace spline
} // end namespace interpolation
} // end namespace numeric

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
numeric::interpolation::spline::SimpleInterpolator::save( Archive & arc ) const {
	arc( cereal::base_class< class numeric::interpolation::spline::Interpolator >( this ) );
	arc( CEREAL_NVP( x_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( y_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( ddy_ ) ); // utility::vector1<Real>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
numeric::interpolation::spline::SimpleInterpolator::load( Archive & arc ) {
	arc( cereal::base_class< class numeric::interpolation::spline::Interpolator >( this ) );
	arc( x_ ); // utility::vector1<Real>
	arc( y_ ); // utility::vector1<Real>
	arc( ddy_ ); // utility::vector1<Real>
}

SAVE_AND_LOAD_SERIALIZABLE( numeric::interpolation::spline::SimpleInterpolator );
CEREAL_REGISTER_TYPE( numeric::interpolation::spline::SimpleInterpolator )

CEREAL_REGISTER_DYNAMIC_INIT( numeric_interpolation_spline_SimpleInterpolator )
#endif // SERIALIZATION
