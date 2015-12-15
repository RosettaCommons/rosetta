// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/interpolation/spline/CompoundInterpolator.cc
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler


#include <numeric/interpolation/spline/CompoundInterpolator.hh>
#include <numeric/interpolation/spline/SimpleInterpolator.hh>
#include <utility/tools/make_vector.hh>

#include <algorithm>

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

struct compare_interp_range {
	bool operator()( interp_range const & a, interp_range const & b ) {
		return a.lb < b.ub;
	}
};

void
CompoundInterpolator::add_range(
	InterpolatorOP interp,
	Real lb,
	Real ub
) {
	interp_range ir;
	ir.lb = lb;
	ir.ub = ub;
	ir.interp = interp;
	interpolators_.push_back( ir );
	std::sort( interpolators_.begin(), interpolators_.end(), compare_interp_range() );
	for ( size_t i = 1; i < interpolators_.size(); ++i ) {
		assert( interpolators_[i].ub <= interpolators_[i+1].lb );
	}
}


void
CompoundInterpolator::interpolate(
	Real x,
	Real & y,
	Real & dy
) {

	if ( has_lb_function() && x < get_lb_function_cutoff() ) {
		return compute_lb_function_solution(x,y);
	}
	if ( has_ub_function() && x > get_ub_function_cutoff() ) {
		return compute_ub_function_solution(x,y);
	}

	for ( size_t i = 1; i <= interpolators_.size(); ++i ) {
		if ( interpolators_[i].lb <= x && x <= interpolators_[i].ub ) {
			return interpolators_[i].interp->interpolate(x,y,dy);
		}
	}
	assert(false);
}

/// @brief serialize the Interpolator to a json_spirit object
utility::json_spirit::Value CompoundInterpolator::serialize()
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;
	std::vector<Value> interpolator_data;
	for ( utility::vector1<interp_range>::iterator it = interpolators_.begin(); it != interpolators_.end(); ++it ) {
		Pair ub("ub",Value(it->ub));
		Pair lb("lb",Value(it->lb));
		Pair interpolator("interp",it->interp->serialize());
		interpolator_data.push_back(Value(utility::tools::make_vector(ub,lb,interpolator)));
	}

	Pair inter_list("interp_list",Value(interpolator_data));
	Pair base_data("base_data",Interpolator::serialize());

	return Value(utility::tools::make_vector(inter_list,base_data));

}

/// @brief deserialize a json_spirit object to a Interpolator
void CompoundInterpolator::deserialize(utility::json_spirit::mObject data)
{
	interpolators_.clear();
	utility::json_spirit::mArray interpolator_data(data["interp_list"].get_array());
	for ( utility::json_spirit::mArray::iterator it = interpolator_data.begin(); it != interpolator_data.end(); ++it ) {
		utility::json_spirit::mObject interpolator_record(it->get_obj());
		InterpolatorOP current_interpolator( new SimpleInterpolator() );
		current_interpolator->deserialize(interpolator_record["interp"].get_obj());
		add_range(current_interpolator,data["lb"].get_real(),data["ub"].get_real());
	}

	Interpolator::deserialize(data["base_data"].get_obj());
}

bool CompoundInterpolator::operator == ( Interpolator const & other ) const
{
	if ( ! Interpolator::operator==( other ) ) return false;

	CompoundInterpolator const & other_downcast( static_cast< CompoundInterpolator const & > ( other ) );
	if ( interpolators_.size() != other_downcast.interpolators_.size() ) return false;
	for ( platform::Size ii = 1; ii <= interpolators_.size(); ++ii ) {
		if ( interpolators_[ii].lb != other_downcast.interpolators_[ii].lb ) return false;
		if ( interpolators_[ii].ub != other_downcast.interpolators_[ii].ub ) return false;
		if ( ! (*interpolators_[ii].interp == *other_downcast.interpolators_[ii].interp) ) return false;
	}
	return true;
}

bool CompoundInterpolator::same_type_as_me( Interpolator const & other ) const
{
	return dynamic_cast< CompoundInterpolator const * > (&other);
}

} // end namespace spline
} // end namespace interpolation
} // end namespace numeric


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
numeric::interpolation::spline::interp_range::save( Archive & arc ) const {
	arc( CEREAL_NVP( lb ) ); // Real
	arc( CEREAL_NVP( ub ) ); // Real
	arc( CEREAL_NVP( interp ) ); // InterpolatorOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
numeric::interpolation::spline::interp_range::load( Archive & arc ) {
	arc( lb ); // Real
	arc( ub ); // Real
	arc( interp ); // InterpolatorOP
}

SAVE_AND_LOAD_SERIALIZABLE( numeric::interpolation::spline::interp_range );

/// @brief Automatically generated serialization method
template< class Archive >
void
numeric::interpolation::spline::CompoundInterpolator::save( Archive & arc ) const {
	arc( cereal::base_class< class numeric::interpolation::spline::Interpolator >( this ) );
	arc( CEREAL_NVP( interpolators_ ) ); // utility::vector1<interp_range>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
numeric::interpolation::spline::CompoundInterpolator::load( Archive & arc ) {
	arc( cereal::base_class< class numeric::interpolation::spline::Interpolator >( this ) );
	arc( interpolators_ ); // utility::vector1<interp_range>
}

SAVE_AND_LOAD_SERIALIZABLE( numeric::interpolation::spline::CompoundInterpolator );
CEREAL_REGISTER_TYPE( numeric::interpolation::spline::CompoundInterpolator )

CEREAL_REGISTER_DYNAMIC_INIT( numeric_interpolation_spline_CompoundInterpolator )
#endif // SERIALIZATION

