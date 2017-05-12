// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/interpolation/Interpolator.hh
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler


#ifndef INCLUDED_numeric_interpolation_spline_Interpolator_hh
#define INCLUDED_numeric_interpolation_spline_Interpolator_hh

#include <numeric/types.hh>

#include <utility/json_spirit/json_spirit_value.h>

#include <utility/pointer/ReferenceCount.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace numeric {
namespace interpolation {
namespace spline {

class Interpolator;
typedef utility::pointer::shared_ptr< Interpolator > InterpolatorOP;
typedef utility::pointer::shared_ptr< Interpolator const > InterpolatorCOP;

class Interpolator : public utility::pointer::ReferenceCount {

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Interpolator();

	Interpolator();

	virtual InterpolatorOP clone() const = 0;

	virtual void interpolate( numeric::Real x, numeric::Real & y, numeric::Real & dy ) const = 0;

	/// @brief set a linear function describing the behavior of the interpolator after a given lower bound.  This lower bound can be distinct from the lb of the spline
	void set_lb_function(numeric::Real const &  lb, numeric::Real const & slope, numeric::Real const & intercept);
	/// @brief set a linear function describing the behavior of the interpolator after a given upper bound.  This upper bound can be distinct from the ub of the spline
	void set_ub_function(numeric::Real const & ub, numeric::Real const & slope,numeric::Real const & intercept);

	/// @brief return true if the interpolator has a defined lower bound function
	bool has_lb_function() const;
	/// @brief return false if the interpolator has a defined upper bound function
	bool has_ub_function() const;

	/// @brief get the lower bound cutoff
	numeric::Real get_lb_function_cutoff() const;
	/// @brief get the upper bound cutoff
	numeric::Real get_ub_function_cutoff() const;

	/// @brief compute the y value of the lower bound function given an x value
	void compute_lb_function_solution(numeric::Real x, numeric::Real & y) const;

	/// @brief compute the y value of the lower bound function given an x value
	void compute_ub_function_solution(numeric::Real x, numeric::Real & y) const;

	/// @brief serialize the Interpolator to a json_spirit object
	virtual utility::json_spirit::Value serialize() const;
	/// @brief deserialize a json_spirit object to a Interpolator
	virtual void deserialize(utility::json_spirit::mObject data);

	virtual bool operator == ( Interpolator const & other ) const;
	virtual bool same_type_as_me( Interpolator const & other ) const;

private:

	bool has_lb_function_;
	bool has_ub_function_;

	numeric::Real lb_cutoff_;
	numeric::Real ub_cutoff_;

	numeric::Real lb_slope_;
	numeric::Real ub_slope_;

	numeric::Real lb_intercept_;
	numeric::Real ub_intercept_;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // end namespace spline
} // end namespace interpolation
} // end namespace numeric

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( numeric_interpolation_spline_Interpolator )
#endif // SERIALIZATION


#endif
