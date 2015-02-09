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

#ifndef INCLUDED_basic_interpolate_hh
#define INCLUDED_basic_interpolate_hh

// Package headers

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray3A.hh>
#include <ObjexxFCL/Fmath.hh>

#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2P.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/StaticIndexRange.fwd.hh>
#include <ObjexxFCL/StaticIndexRange.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <iosfwd>
#include <limits>
#include <string>

// util_interpolate Function Declarations
namespace basic {

void
interpolate_get_angle_bins(
	double const x,
	double const binrange,
	int & xbin,
	int & xbin_next,
	double & xd
);


void
interpolate_bilinear(
	int const xbin,
	int const xbin_next,
	double const xd,
	int const ybin,
	int const ybin_next,
	double const yd,
	ObjexxFCL::FArray2A< double > xy_func,
	int const xbin_count,
	int const ybin_count,
	double const binrange,
	bool const angles,
	double & val,
	double & dval_dx,
	double & dval_dy
);


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
);

void
interpolate_trilinear(
	int const xbin,
	int const xbin_next,
	double const xd,
	int const ybin,
	int const ybin_next,
	double const yd,
	int const zbin,
	int const zbin_next,
	double const zd,
	ObjexxFCL::FArray3A< double > xyz_func,
	int const xbin_count,
	int const ybin_count,
	int const zbin_count,
	double const binrange, // assumes that have the same bin size in each dimention
	bool const angles,
	double & val,
	double & dval_dx,
	double & dval_dy,
	double & dval_dz
);

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
	double const binrange, // assumes that have the same bin size in each dimention
	bool const angles,
	double & val,
	double & dval_dx,
	double & dval_dy,
	double & dval_dz
);

void
interpolate_2d_func_of_angles(
	double const x,
	double const y,
	ObjexxFCL::FArray2A< double > xy_func,
	double & val,
	double & dval_dx,
	double & dval_dy
);


////////////////////////////////////////////////////////////////////////////////
/// @begin interpolate_get_angle_bins
///
/// @brief get bin information for a periodic value w/ periodic bins
///
/// @detailed
///
///     for 'x', an angle in degrees, and angular bins of width
///     'binrange' aligned to start bin 1 at a value of 0. degrees, find
///     the two bins, whose average bin values 'x' falls between, and
///     report the error between the lower average bin value and the real
///     value of 'x'
///
/// @param[in]   x - in - angle in degrees
/// @param[in]   binrange - in - degrees per bin
/// @param[out]   xbin - out - anglular bin whose average value is the lower
///                      of the two average bin values bracketing 'x'
/// @param[out]   xbin_next - out - anglular bin whose average value is the lower
///                      of the two average bin values bracketing 'x'
/// @param[out]   xd - out - error between the average value of bin 'xbin' and 'x'
///
/// @remarks
///
/// @references
///
/// @authors
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
inline
void
interpolate_get_angle_bins(
	double const x,
	double const binrange,
	int & xbin,
	int & xbin_next,
	double & xd
)
{
//------------------------------------------------------------------------------
// ctsa - get nbins for binrange
//
	int const nbins = static_cast< int >( 360.0 / binrange );

// ctsa -  convert angle to double(anglebin) and
//         mod value to range: [1.,double(nbins)+1.)
//
//   note that:  double(anglebin) = 1. + (angle-.5*binrange)/binrange
//     ...this solves for the lowest of the two bracketing bin averages,
//       rather  than the nearest bin average, and is thus only
//       appropriate for interpolation
//

	double const xbin_real = 1.0 + ObjexxFCL::mod(
	 ObjexxFCL::mod( (x-(0.5*binrange))/binrange, static_cast< double >(nbins) ) + nbins,
	 static_cast< double >(nbins) );

// ctsa -  convert double bin values to array lookup bin values
//
	xbin = static_cast< int >( xbin_real );

// ctsa -  get next lookup bin and convert to range: (1,nbins)
//
	xbin_next = 1 + ObjexxFCL::mod( xbin, nbins );

// ctsa -  get error
//
	xd = xbin_real - xbin;


	//// ctsa -- 5-2003
	////   ...another double resolution bug of type:
	////
	//// static_cast< int >(double(intval)-epsilon) == intval, epsilon << 1
	////
	//// this causes xbin to sporatically resolve to bins+1, a
	//// very bad thing; check below should fix...
	////
	if ( xbin == nbins + 1 ) {
		xbin = 1;
		xd = 0.0;
	}
}


////////////////////////////////////////////////////////////////////////////////
/// @begin interpolate_bilinear
///
/// @brief get bilinear interpolate of a 2d periodic function
///
/// @detailed
///
///     Value and derivatives of bilinear interpolation on a 2d binned
///     periodic function, represented as a 2d array.  (see Num. Recipes
///     v2, sec 3.6, "Interpolation in two or more dimensions")
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
/// @param[in]   xbin - in - low bin in 1st dim
/// @param[in]   xbin_next - in - high bin in 1st dim
/// @param[in]   xd - in - error between low and real val in 1st dim
/// @param[in]   ybin - in - " 2nd dim "
/// @param[in]   ybin_next - in - " 2nd dim "
/// @param[in]   yd - in - " 2nd dim "
/// @param[in]   xy_func - in - 2d-array specifying the 2d binned function values,
///                        assumed to be periodic in each dimension
/// @param[in]   xbin_count - in - number of bins in the 1st dim of the periodic
///                        function
/// @param[in]   ybin_count - in - " 2nd dim "
/// @param[in]   binrange - in - range of bin in angles ( for both dimensions )
/// @param[in]   angles - in - if true, treat the the *values* of xy_func as
///                       as having a periodicity of 360. ( note that
///                       the bin ranges are already assumed periodic )
/// @param[out]   val - out - bilinear interpolated value of xy_func
/// @param[out]   dval_dx - out - derivative of interpolation w.r.t 1st dim
/// @param[out]   dval_dy - out - " 2nd dim "
///
/// @remarks
///
/// @references
///
/// @authors ctsa 8-19-03
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
inline
void
interpolate_bilinear(
	int const xbin,
	int const xbin_next,
	double const xd,
	int const ybin,
	int const ybin_next,
	double const yd,
	ObjexxFCL::FArray2A< double > xy_func,
	int const xbin_count,
	int const ybin_count,
	double const binrange,
	bool const angles,
	double & val,
	double & dval_dx,
	double & dval_dy
)
{
	xy_func.dimension( xbin_count, ybin_count );

//------------------------------------------------------------------------------
	double const x0y0 = xy_func(xbin,ybin);
	double const x1y0 = xy_func(xbin_next,ybin);
	double const x0y1 = xy_func(xbin,ybin_next);
	double const x1y1 = xy_func(xbin_next,ybin_next);

	interpolate_bilinear_by_value(x0y0,x1y0,x0y1,x1y1,xd,yd,binrange,angles,val,
	 dval_dx,dval_dy);
}

inline
void
interpolate_trilinear(
	int const xbin,
	int const xbin_next,
	double const xd,
	int const ybin,
	int const ybin_next,
	double const yd,
	int const zbin,
	int const zbin_next,
	double const zd,
	ObjexxFCL::FArray3A< double > xyz_func,
	int const xbin_count,
	int const ybin_count,
	int const zbin_count,
	double const binrange, // assumes that have the same bin size in each dimention
	bool const angles,
	double & val,
	double & dval_dx,
	double & dval_dy,
	double & dval_dz
)
{
	xyz_func.dimension( xbin_count, ybin_count, zbin_count );

	double const x0y0z0 = xyz_func( xbin, ybin, zbin );
	double const x1y0z0 = xyz_func( xbin_next, ybin, zbin );
	double const x0y1z0 = xyz_func( xbin, ybin_next, zbin );
	double const x1y1z0 = xyz_func( xbin_next,ybin_next, zbin );
	double const x0y0z1 = xyz_func( xbin, ybin, zbin_next );
	double const x1y0z1 = xyz_func( xbin_next, ybin, zbin_next );
	double const x0y1z1 = xyz_func( xbin, ybin_next, zbin_next );
	double const x1y1z1 = xyz_func( xbin_next, ybin_next, zbin_next );

	interpolate_trilinear_by_value(
		x0y0z0, x1y0z0, x0y1z0, x1y1z0,
		x0y0z1, x1y0z1, x0y1z1, x1y1z1,
		xd, yd, zd, binrange, angles,
		val, dval_dx, dval_dy, dval_dz );
}

////////////////////////////////////////////////////////////////////////////////
/// @begin interpolate_2d_func_of_angles
///
/// @brief get bilinear interpolate of a 2d function with degree angle arguments
///
/// @detailed
///
///     Value and derivatives of bilinear interpolation on a 2d function
///     with degree angle arguments represented by a 2d array with binned
///     degree angle indices.  (see Num. Recipes v2, sec 3.6,
///     "Interpolation in two or more dimensions")
///
///     The derivatives are of the bilinear interpolation only. This is
///     not the same value as a numerical derivative of function
///     values. The derivatives here are discontinuous every time you
///     cross over a bin boundary in either dimension
///
/// @param[in]   x - in - function argument value in 1st dim to find interpolate for
/// @param[in]   y - in - " 2nd dim "
/// @param[in]   xy_func - in - 2d array with binned angle value indicies;
///               represents known values of 2d function being interpolated.
/// @param[out]   val - out - bilinear interpolated value
/// @param[out]   dval_dx - out - derivative of interpolated value w.r.t. 1st dim angle
/// @param[out]   dval_dy - out - " 2nd dim "
///
/// @remarks
///
///     !!! assumes 10 degree bins in both dimensions of xy_func
///
/// @references
///
/// @authors ctsa 8-19-03
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
inline
void
interpolate_2d_func_of_angles(
	double const x,
	double const y,
	ObjexxFCL::FArray2A< double > xy_func,
	double & val,
	double & dval_dx,
	double & dval_dy
)
{
	int const nbins = { 36 };
	xy_func.dimension( nbins, nbins );

//ctsa fixed parameters
	double const binrange = { 10.0f };

//ctsa local
	int xbin,ybin,xbin_next,ybin_next;
	double xd,yd;

//------------------------------------------------------------------------------
	interpolate_get_angle_bins(x,binrange,xbin,xbin_next,xd);
	interpolate_get_angle_bins(y,binrange,ybin,ybin_next,yd);

	bool const treat_as_angles = false;

	interpolate_bilinear(xbin,xbin_next,xd,ybin,ybin_next,yd,xy_func,nbins,nbins,
	 binrange,treat_as_angles,val,dval_dx,dval_dy);
}


} // basic

#endif
