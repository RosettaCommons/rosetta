// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/interpolation/Interpolator.hh
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler
///

#ifndef INCLUDED_numeric_interpolation_spline_Interpolator_hh
#define INCLUDED_numeric_interpolation_spline_Interpolator_hh

#include <numeric/types.hh>

#include <utility/vector1.hh>
#include <utility/json_spirit/json_spirit_value.h>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>

namespace numeric {
namespace interpolation {
namespace spline {

using numeric::Real;

class Interpolator : public utility::pointer::ReferenceCount {

public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~Interpolator();

	Interpolator();

	virtual void interpolate( Real x, Real & y, Real & dy ) = 0;

	/// @brief set a linear function describing the behavior of the interpolator after a given lower bound.  This lower bound can be distinct from the lb of the spline
	void set_lb_function(Real const &  lb, Real const & slope, Real const & intercept);
	/// @brief set a linear function describing the behavior of the interpolator after a given upper bound.  This upper bound can be distinct from the ub of the spline
	void set_ub_function(Real const & ub, Real const & slope,Real const & intercept);

	/// @brief return true if the interpolator has a defined lower bound function
	bool has_lb_function() const;
	/// @brief return false if the interpolator has a defined upper bound function
	bool has_ub_function() const;

	/// @brief get the lower bound cutoff
	Real get_lb_function_cutoff() const;
	/// @brief get the upper bound cutoff
	Real get_ub_function_cutoff() const;

	/// @brief compute the y value of the lower bound function given an x value
	void compute_lb_function_solution(Real x, Real & y) const;

	/// @brief compute the y value of the lower bound function given an x value
	void compute_ub_function_solution(Real x, Real & y) const;

	/// @brief serialize the Interpolator to a json_spirit object
	virtual utility::json_spirit::Value serialize();
	/// @brief deserialize a json_spirit object to a Interpolator
	virtual void deserialize(utility::json_spirit::mObject data);

private:

	bool has_lb_function_;
	bool has_ub_function_;

	Real lb_cutoff_;
	Real ub_cutoff_;

	Real lb_slope_;
	Real ub_slope_;

	Real lb_intercept_;
	Real ub_intercept_;


};

typedef utility::pointer::shared_ptr< Interpolator > InterpolatorOP;

} // end namespace spline
} // end namespace interpolation
} // end namespace numeric

#endif
