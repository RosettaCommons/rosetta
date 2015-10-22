// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author


// Rosetta Headers
#include <basic/interpolate.hh>
#include <basic/basic.hh>  // for subtract_degree_angles, angle_in_range

namespace basic {

////////////////////////////////////////////////////////////////////////////////
///
/// @brief get bilinear interpolate, given four points of a 2d periodic function
///
/// @details
///
///     Value and derivatives of bilinear interpolation between four
///     specified input values, which are the function values
///     corresponding to four points which bracket the interpolated point
///     in 2-dimensions.  (see diagram below) (see Num. Recipes v2, sec
///     3.6, "Interpolation in two or more dimensions")
///
///     Note that it is automatically assumed that the arguments are
///     periodic, that is f(a,b),a,b=periodic value; the *value* of the
///     function can be treated as periodic ( fixed to a periodicity of
///     360.) as well when the input parameter angle is set to true
///
///     The derivatives are of the bilinear interpolation only. This is
///     not the same value as a numerical derivative of function
///     values. The derivatives here are discontinuous every time you
///     cross over a bin boundary in either dimension
///
///
/// -  bilinear interpolation:
///
///    x0y1 +----------+ x1y1
///         |          |
///   1-yd  |    (x,y) |
///         - - +      |
///     yd  |   .      |
///         |   .      |
///    x0y0 +---|------+ x1y0
///          xd   1-xd
///
///
///
/// @param[in]   x0y0 - in - input bracketing function value (see diagram)
/// @param[in]   x1y0 - in - " "
/// @param[in]   x0y1 - in - " "
/// @param[in]   x1y1 - in - " "
/// @param[in]   xd - in - error value in the 1st dim between lower lst dim
///                   bracket value 'x0' and point to be interpolated, 'x'
/// @param[in]   yd - in - " 2nd dim "
/// @param[in]   binrange - in - range of bin in angles ( for both dimensions )
/// @param[in]   angles - in - if true, treat the the *values* of xy_func as
///                       as having a periodicity of 360. ( note that
///                       the bin ranges are already assumed periodic )
/// @param[out]   val - out - biliinear interpolated value
/// @param[out]   dval_dx - out - derivative of 'val' in the 1st ('x') dim
/// @param[out]   dval_dy - out - " 2nd dim "
///
/// @remarks
///
/// @references
///
/// @author ctsa 8-19-03
///
/////////////////////////////////////////////////////////////////////////////////
void
interpolate_bilinear_by_value(
	double const x0y0,
	double const x1y0,
	double const x0y1,
	double const x1y1,
	double const xd,
	double const yd,
	double const binrange,
	bool const angles,
	double & val,
	double & dval_dx,
	double & dval_dy
)
{
	double const w1 = ( 1.0f - xd ) * ( 1.0f - yd );
	double const w2 = xd * ( 1.0f - yd );
	double const w3 = ( 1.0f - xd ) * yd;
	double const w4 = xd * yd;

	if ( angles ) {
		//// ctsa - hack method for weighted angle averaging, if we could burn
		////  cpu time, then the right thing to do would probably be to
		////  break func into x and y components, add up, interpolate,
		////  and find the direction

		double const w12 = w1 + w2;
		double a12;
		if ( w12 != 0.0 ) {
			a12 = ( w1 * x0y0 + w2 * ( subtract_degree_angles(x1y0,x0y0) + x0y0 ) ) / w12;
		} else {
			a12 = 0.0;
		}

		double const w34 = w3 + w4;
		double a34;
		if ( w34 != 0.0 ) {
			a34 = ( w3 * x0y1 + w4 * ( subtract_degree_angles(x1y1,x0y1) + x0y1 ) ) / w34;
		} else {
			a34 = 0.0;
		}

		val = w12 * a12 + w34 * ( subtract_degree_angles(a34,a12) + a12 );
		angle_in_range(val);

		dval_dx = ( ( 1.0f - yd ) * subtract_degree_angles(x1y0,x0y0) + yd * subtract_degree_angles(x1y1,x0y1) ) / binrange;
		dval_dy = ( ( 1.0f - xd ) * subtract_degree_angles(x0y1,x0y0) + xd * subtract_degree_angles(x1y1,x1y0) ) / binrange;
	} else {
		val = x0y0 * w1 + x1y0 * w2 + x0y1 * w3 + x1y1 * w4;

		dval_dx = ( ( 1.0f - yd ) * ( x1y0 - x0y0 ) + yd * ( x1y1 - x0y1 ) ) / binrange;
		dval_dy = ( ( 1.0f - xd ) * ( x0y1 - x0y0 ) + xd * ( x1y1 - x1y0 ) ) / binrange;
	}
}

/// @brief paraphazed from above for three dimentions
/// Without additional scaling this will only work if
/// the bin sizes in each dimention are equal I think.
/// Otherwise it could be computed as the product of
/// three linear interpolations
///
/// Perform linear interpolation between x0y0z0 and x1y0z0 to find y0z0
/// Preform linear interpolation between x0y0z1 and x1y0z1 to find y0z1
/// Preform linear interpolation between x0y1z1 and x1y1z1 to find y1z1
/// Preform linear interpolation between x0y1z0 and x1y1z0 to find y1z0
/// Preform linear interpolation between y0z0 and y1z0 to find z0
/// Preform linear interpolation between y0z1 and y1z1 to find z1
/// Preform linear interpolation between z0 and z1 to find val

void
interpolate_trilinear_by_value(
	double const x0y0z0,
	double const x1y0z0,
	double const x0y1z0,
	double const x1y1z0,
	double const x0y0z1,
	double const x1y0z1,
	double const x0y1z1,
	double const x1y1z1,
	double const xd,
	double const yd,
	double const zd,
	double const binrange,
	bool const angles,
	double & val,
	double & dval_dx,
	double & dval_dy,
	double & dval_dz
)
{
	double const w1 = ( 1.0f - xd ) * ( 1.0f - yd ) * ( 1.0f -zd );
	double const w2 = xd * ( 1.0f - yd ) * ( 1.0f -zd );
	double const w3 = ( 1.0f - xd ) * yd * ( 1.0f -zd );
	double const w4 = ( 1.0f - xd ) * ( 1.0f - yd ) * zd;
	double const w5 = xd * yd * ( 1.0f -zd );
	double const w6 = xd * ( 1.0f - yd ) * zd;
	double const w7 = ( 1.0f - xd ) * yd * zd;
	double const w8 = xd * yd * zd;

	if ( angles ) {
		// not sure why angles need to be treated differently
		val = x0y0z0 * w1 + x1y0z0 * w2 + x0y1z0 * w3 + x0y0z1 * w4 + x1y1z0 * w5 + x1y0z1 * w6 + x0y1z1 * w7 + x1y1z1 * w8;
		dval_dx = ( -1.0f * x0y0z0 * ( 1.0f - yd ) * ( 1.0f - zd ) + x1y0z0 * ( 1.0f - yd ) * ( 1.0f - zd ) - x0y1z0 * yd * ( 1.0f - zd ) + x1y0z1 * yd * ( 1.0f - zd ) - x1y0z1 * ( 1.0f - yd ) * zd + x1y0z0 * ( 1.0f - yd ) * zd - x1y0z0 * yd * zd + x1y1z1 * yd * zd  ) / binrange;
		dval_dy = ( -1.0f * x0y0z0 * ( 1.0f - xd ) * ( 1.0f - zd ) + x0y1z0 * ( 1.0f - xd ) * ( 1.0f - zd ) - x1y0z0 * xd * ( 1.0f - zd ) + x1y0z1 * xd * ( 1.0f - zd ) - x1y0z1 * ( 1.0f - xd ) * zd + x1y0z0 * ( 1.0f - xd ) * zd - x1y0z0 * xd * zd + x1y1z1 * xd * zd  ) / binrange;
		dval_dz = ( -1.0f * x0y0z0 * ( 1.0f - xd ) * ( 1.0f - yd ) + x1y0z1 * ( 1.0f - xd ) * ( 1.0f - yd ) - x1y0z0 * xd * ( 1.0f - yd ) + x1y0z0 * xd * ( 1.0f - yd ) - x0y1z0 * ( 1.0f - xd ) * yd + x1y0z0 * ( 1.0f - xd ) * yd - x1y0z1 * xd * yd + x1y1z1 * xd * yd  ) / binrange;

	} else {
		val = x0y0z0 * w1 + x1y0z0 * w2 + x0y1z0 * w3 + x0y0z1 * w4 + x1y0z1 * w5 + x0y1z1 * w6 + x1y1z0 * w7 + x1y1z1 * w8;
		dval_dx = ( -1.0f * x0y0z0 * ( 1.0f - yd ) * ( 1.0f - zd ) + x1y0z0 * ( 1.0f - yd ) * ( 1.0f - zd ) - x0y1z0 * yd * ( 1.0f - zd ) + x1y0z1 * yd * ( 1.0f - zd ) - x1y0z1 * ( 1.0f - yd ) * zd + x1y0z0 * ( 1.0f - yd ) * zd - x1y0z0 * yd * zd + x1y1z1 * yd * zd  ) / binrange;
		dval_dy = ( -1.0f * x0y0z0 * ( 1.0f - xd ) * ( 1.0f - zd ) + x0y1z0 * ( 1.0f - xd ) * ( 1.0f - zd ) - x1y0z0 * xd * ( 1.0f - zd ) + x1y0z1 * xd * ( 1.0f - zd ) - x1y0z1 * ( 1.0f - xd ) * zd + x1y0z0 * ( 1.0f - xd ) * zd - x1y0z0 * xd * zd + x1y1z1 * xd * zd  ) / binrange;
		dval_dz = ( -1.0f * x0y0z0 * ( 1.0f - xd ) * ( 1.0f - yd ) + x1y0z1 * ( 1.0f - xd ) * ( 1.0f - yd ) - x1y0z0 * xd * ( 1.0f - yd ) + x1y0z0 * xd * ( 1.0f - yd ) - x0y1z0 * ( 1.0f - xd ) * yd + x1y0z0 * ( 1.0f - xd ) * yd - x1y0z1 * xd * yd + x1y1z1 * xd * yd  ) / binrange;

	}
}

} // namespace basic
