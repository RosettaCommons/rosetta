// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/interpolation/spline/SimpleInterpolator.hh
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler
///

#ifndef INCLUDED_numeric_interpolation_spline_SimpleInterpolator_hh
#define INCLUDED_numeric_interpolation_spline_SimpleInterpolator_hh

#include <numeric/interpolation/spline/Interpolator.hh>
#include <utility/json_spirit/json_spirit_writer.h>
#include <utility/json_spirit/json_spirit_reader.h>

#include <numeric/interpolation/spline/spline_functions.hh>

namespace numeric {
namespace interpolation {
namespace spline {

using numeric::Real;
using utility::vector1;


class SimpleInterpolator : public Interpolator {

public:

	SimpleInterpolator(
	  utility::vector1<Real> const & x,
	  utility::vector1<Real> const & y,
	  Real lbdy,
	  Real ubdy
	);

	SimpleInterpolator();

	void interpolate( Real x, Real & y, Real & dy );

	/// @brief serialize the Interpolator to a json_spirit object
	virtual utility::json_spirit::Value serialize();
	/// @brief deserialize a json_spirit object to a Interpolator
	virtual void deserialize(utility::json_spirit::mObject data);

private:

	utility::vector1<Real> x_, y_, ddy_;

};

typedef utility::pointer::owning_ptr< Interpolator > InterpolatorOP;

} // end namespace spline
} // end namespace interpolation
} // end namespace numeric

#endif
