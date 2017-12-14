// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/interpolation/spline/Interpolator.cc
/// @author Sam DeLuca

#include <numeric/interpolation/spline/Interpolator.hh>

#include <utility/tools/make_vector.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace numeric {
namespace interpolation {
namespace spline {

/// @details Auto-generated virtual destructor
Interpolator::~Interpolator() = default;

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


utility::json_spirit::Value Interpolator::serialize() const
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair lb_fxn("lbfxn",Value(has_lb_function_));
	Pair ub_fxn("ubfxn",Value(has_ub_function_));

	Pair lb_cut("lbcut",Value(lb_cutoff_));
	Pair ub_cut("ubcut",Value(ub_cutoff_));

	Pair lb_slope("lbslope",Value(lb_slope_));
	Pair ub_slope("ubslope",Value(ub_slope_));

	Pair lb_int("lbint",Value(lb_intercept_));
	Pair ub_int("ubint",Value(ub_intercept_));

	return Value(utility::tools::make_vector(lb_fxn,ub_fxn,lb_cut,ub_cut,lb_slope,ub_slope,lb_int,ub_int));
}

void Interpolator::deserialize(utility::json_spirit::mObject data)
{
	has_lb_function_ = data["lbfxn"].get_bool();
	has_ub_function_ = data["ubfxn"].get_bool();

	lb_cutoff_ = data["lbcut"].get_real();
	ub_cutoff_ = data["ubcut"].get_real();

	lb_slope_ = data["lbslope"].get_real();
	ub_slope_ = data["ubslope"].get_real();

	lb_intercept_ = data["lbint"].get_real();
	ub_intercept_ = data["ubint"].get_real();

}

bool Interpolator::operator == ( Interpolator const & other ) const
{
	if ( !       same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	if ( has_lb_function_   != other.has_lb_function_   ) return false;
	if ( has_ub_function_   != other.has_ub_function_   ) return false;
	if ( lb_cutoff_         != other.lb_cutoff_         ) return false;
	if ( ub_cutoff_         != other.ub_cutoff_         ) return false;
	if ( lb_slope_          != other.lb_slope_          ) return false;
	if ( ub_slope_          != other.ub_slope_          ) return false;
	if ( lb_intercept_      != other.lb_intercept_      ) return false;
	if ( ub_intercept_      != other.ub_intercept_      ) return false;

	return true;
}

bool Interpolator::same_type_as_me( Interpolator const & ) const
{
	return true;
}

}
}
}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
numeric::interpolation::spline::Interpolator::save( Archive & arc ) const {
	arc( CEREAL_NVP( has_lb_function_ ) ); // _Bool
	arc( CEREAL_NVP( has_ub_function_ ) ); // _Bool
	arc( CEREAL_NVP( lb_cutoff_ ) ); // Real
	arc( CEREAL_NVP( ub_cutoff_ ) ); // Real
	arc( CEREAL_NVP( lb_slope_ ) ); // Real
	arc( CEREAL_NVP( ub_slope_ ) ); // Real
	arc( CEREAL_NVP( lb_intercept_ ) ); // Real
	arc( CEREAL_NVP( ub_intercept_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
numeric::interpolation::spline::Interpolator::load( Archive & arc ) {
	arc( has_lb_function_ ); // _Bool
	arc( has_ub_function_ ); // _Bool
	arc( lb_cutoff_ ); // Real
	arc( ub_cutoff_ ); // Real
	arc( lb_slope_ ); // Real
	arc( ub_slope_ ); // Real
	arc( lb_intercept_ ); // Real
	arc( ub_intercept_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( numeric::interpolation::spline::Interpolator );
CEREAL_REGISTER_TYPE( numeric::interpolation::spline::Interpolator )

CEREAL_REGISTER_DYNAMIC_INIT( numeric_interpolation_spline_Interpolator )
#endif // SERIALIZATION
