// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/numeric/interpolation/SplineGenerator.hh
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler
///

#ifndef INCLUDED_numeric_interpolation_spline_SplineGenerator_hh
#define INCLUDED_numeric_interpolation_spline_SplineGenerator_hh

#include <numeric/interpolation/spline/Interpolator.hh>

#include <numeric/types.hh>
#include <utility/vector1.hh>

namespace numeric {
namespace interpolation {
namespace spline {

using numeric::Real;
using utility::vector1;

struct Point {
	Point( Real xin, Real yin            ) : x(xin), y(yin), dy(-12345.0), has_dy(false) {}
	Point( Real xin, Real yin, Real dyin ) : x(xin), y(yin), dy(  dyin  ), has_dy(true ) {}
	Real x;
	Real y;
	Real dy;
	bool has_dy;
};

class SplineGenerator {
public:

	SplineGenerator(
		Real lbx, Real lby, Real lbdy,
		Real ubx, Real uby, Real ubdy
	);

	SplineGenerator();

	~SplineGenerator();
	void add_known_value( Real x, Real y );

	void add_known_value( Real x, Real y, Real dy );

	InterpolatorOP get_interpolator();

	//getters for the rest of the private data
	Real get_lbx()
	{
		return lbx_;
	}

	Real get_lby()
	{
		return lby_;
	}

	Real get_lbdy()
	{
		return lbdy_;
	}

	Real get_ubx()
	{
		return ubx_;
	}

	Real get_uby()
	{
		return uby_;
	}

	Real get_ubdy()
	{
		return ubdy_;
	}

	numeric::Size get_num_points()
	{
		return points_.size();
	}

private:

	Real lbx_,lby_,lbdy_,ubx_,uby_,ubdy_;

	vector1<Point> points_;

	InterpolatorOP interpolator_;

};

} // end namespace spline
} // end namespace interpolation
} // end namespace numeric

#endif
