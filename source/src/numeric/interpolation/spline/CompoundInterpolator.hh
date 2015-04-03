// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/interpolation/spline/CompoundInterpolator.hh
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler


#ifndef INCLUDED_numeric_interpolation_spline_CompoundInterpolator_hh
#define INCLUDED_numeric_interpolation_spline_CompoundInterpolator_hh

#include <numeric/types.hh>

#include <numeric/interpolation/spline/Interpolator.hh>

#include <utility/vector1.hh>

namespace numeric {
namespace interpolation {
namespace spline {

class CompoundInterpolator;
typedef utility::pointer::shared_ptr< CompoundInterpolator > CompoundInterpolatorOP;
typedef utility::pointer::shared_ptr< CompoundInterpolator const > CompoundInterpolatorCOP;

using numeric::Real;

struct interp_range {
	Real lb;
	Real ub;
	InterpolatorOP interp;
};

class CompoundInterpolator : public Interpolator {

public:

	void add_range( InterpolatorOP interp, Real lb, Real ub );

	void interpolate( Real x, Real & y, Real & dy );

	/// @brief serialize the Interpolator to a json_spirit object
	virtual utility::json_spirit::Value serialize();
	/// @brief deserialize a json_spirit object to a Interpolator
	virtual void deserialize(utility::json_spirit::mObject data);

private:

	utility::vector1< interp_range > interpolators_;

};

} // end namespace spline
} // end namespace interpolation
} // end namespace numeric

#endif
