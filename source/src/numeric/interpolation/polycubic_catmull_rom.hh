// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/interpolation/polycubic_catmull_rom.cc
/// @author Rhiju Das

#ifndef INCLUDED_numeric_interpolation_polycubic_catmull_rom
#define INCLUDED_numeric_interpolation_polycubic_catmull_rom

#include <utility/fixedsizearray1.hh>
#include <utility/vector1.hh>
#include <utility/modulo.hh>

#include <numeric/MathNTensor.hh>
#include <cmath>

///////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details
///
///   Also called Keys spline or polycubic 'convolution'.
///   See: https://en.wikipedia.org/wiki/Bicubic_interpolation
///
///   Unlike prior Rosetta spline interpolation (see interpolation::polycubic_interpolation
///     or B-spline implementation in core::scoring::electron_density::SplineInterp),
///     this does *not* require pre-training of
///     spline which can be both time-intensive & memory intensive. Instead, uses
///     an interpolant based on nearest neighbor grid points (like polylinear
///     interpolation) and next-nearest neighbor.
///
///   First derivative is continuous (but second deriv is not, in general).
///
///   If we need more accuracy & smoothness, there is a higher-order cubic spline that
///     looks out one more neighbor (see Keys, "Cubic convolution interpolation
///     for digital image processing",  IEEE Transactions on Acoustics, Speech,
///     and Signal Processing 1981).
///
///  If we need more speed, following is pretty inefficient to ensure generality (see
///     notes). Also, there is apparently a faster version called osculatory interpolation
///     that predated Keys & Catmull & Rom by 80 years; see Meijering & Unser,
///     "A Note on Cubic Convolution Interpolation", IEEE Transactions on Image Processing 2003.
///
///////////////////////////////////////////////////////////////////////////////////////////////////////

namespace numeric {
namespace interpolation {

// @brief for catmull_rom cubic interpolation: boundary condition options for each dimension of tensor
enum CatmullRomSplineBoundaryType
{
	PERIODIC, // define grid points on -180, -170,.. 170, and next value (180) is assumed to be -180.
	FLAT,     // out-of-bounds points are mapped to nearest defined neighbor
	LINEAR    // out-of-bounds points are linearly extrapolated
	// CUBIC // not implemented
};

// @brief enable text output of enum CatmullRomSplineBoundaryType
inline
std::string
to_string( CatmullRomSplineBoundaryType const & type )
{
	switch ( type ) {
	case PERIODIC :
		return "PERIODIC";
	case FLAT :
		return "FLAT";
	case LINEAR :
		return "LINEAR";
	}
	return "UNKNOWN";
}

// @brief  Catmull-Rom interpolation in 1 dimension
template< typename T >
inline
T
catmull_rom_interpolate_basic( utility::fixedsizearray1< T, 4 > const & p,
	Real const & x,
	Real & dval_dx,
	utility::fixedsizearray1< T, 4 > & dval_dp,
	bool const compute_deriv = true )
{
	T const x2 = x*x;
	T const x3 = x2*x;
	T const val =
		(   -0.5*x +     x2 - 0.5*x3)  * p[1] +
		( 1        - 2.5*x2 + 1.5*x3 ) * p[2] +
		(     0.5*x +  2*x2 - 1.5*x3)  * p[3] +
		(           -0.5*x2 + 0.5*x3)  * p[4];

	if ( compute_deriv ) {
		dval_dx =
			(   -0.5 +    2*x - 1.5*x2)  * p[1] +
			(          -  5*x + 4.5*x2 ) * p[2] +
			(     0.5 +   4*x - 4.5*x2)  * p[3] +
			(            -1*x + 1.5*x2)  * p[4];

		dval_dp[1] =    -0.5*x +     x2 - 0.5*x3;
		dval_dp[2] =  1        - 2.5*x2 + 1.5*x3;
		dval_dp[3] =      0.5*x +  2*x2 - 1.5*x3;
		dval_dp[4] =            -0.5*x2 + 0.5*x3;
	}

	return val;
}

// @brief  Catmull-Rom interpolation in 1 dimension
template< typename T >
inline
T
catmull_rom_interpolate_basic( utility::fixedsizearray1< T, 4 > const & p,
	T const & x,
	T & dval_dx,
	bool const compute_deriv = true )
{
	utility::fixedsizearray1< T, 4 > dval_dp;
	return catmull_rom_interpolate_basic( p, x, dval_dx, dval_dp, compute_deriv );
}

// @brief  Catmull-Rom interpolation in 1 dimension
template< typename T >
inline
T
catmull_rom_interpolate_basic( utility::fixedsizearray1< T, 4 >  const & p,
	T const & x )
{
	T dval_dx;
	utility::fixedsizearray1< T, 4 > dval_dp;
	return catmull_rom_interpolate_basic( p, x, dval_dx, false );
}

// @brief  Catmull-Rom interpolation in N dimensions
// @details
//   Interpolates into the tensor F_patch at point d,
//   using a recursion on the number of dimensions, ending
//   in the 1D interpolation catmull_rom_interpolate_basic
//
// TODO: send around F_patch as a const & MathNTensor and pass which_dim to
//        recursively go down each dimension (should save lots of time).
//
// @author rhiju
template< typename T >
inline
T
catmull_rom_interpolate( utility::vector1< T > const & F_patch,
	utility::vector1< T > const & d,
	utility::vector1< T > & deriv,
	bool const & compute_deriv )
{
	T val( 0.0 ), dv_dx( 0.0 );
	utility::fixedsizearray1< T, 4 > p( 0.0 );
	Size const N( d.size() );
	if ( N == 1 ) {
		// TODO this copies code with next chunk below -- unify.
		for ( Size n = 1; n <= 4; n++ ) p[ n ] = F_patch[ n ];
		val = catmull_rom_interpolate_basic( p, d[ 1 ], dv_dx, compute_deriv );
		deriv[ 1 ] = dv_dx;
	} else {
		utility::vector1< T > d_next( N - 1 );
		for ( Size n = 1; n <= N - 1; n++ ) d_next[ n ] = d[ n + 1 ];
		utility::fixedsizearray1< utility::vector1< T >, 4 > dp_dy( utility::vector1< T >( N-1, 0.0 ) );
		Size next_patch_size( pow( 4, N-1 ) );
		utility::vector1< T > F_subpatch( next_patch_size, 0.0 );
		for ( Size i = 1; i <= 4; i++ ) {
			// need to look at i-th sub-tensor of F_patch with 1 fewer dimension
			//  Take advantage of the way the numbers are 'laid out' in F_patch.
			//
			// Would be faster to *not* do an explicit copy but to instead provide a 'vector1' that
			//  points to the same addresses in memory.
			for ( Size n = 1; n <= next_patch_size; n++ ) F_subpatch[ n ] = F_patch[ n + (i-1) * next_patch_size ];
			p[ i ] = catmull_rom_interpolate( F_subpatch, d_next, dp_dy[ i ], compute_deriv );
		}
		utility::fixedsizearray1< T, 4 > dv_dp( 0.0 );
		utility::vector1< T > dv_dy( N-1, 0.0 );
		val = catmull_rom_interpolate_basic( p, d[ 1 ], dv_dx, dv_dp, compute_deriv );
		if ( compute_deriv ) {
			for ( Size n = 1; n <= N-1; n++ ) {
				for ( Size i = 1; i <= 4; i++ ) {
					dv_dy[ n ] += dv_dp[ i ] * dp_dy[ i ][ n ]; // chain rule
				}
			}
		}
		deriv[ 1 ] = dv_dx;
		if ( compute_deriv ) {
			for ( Size n = 1; n <= N-1; n++ ) deriv[ n + 1 ] = dv_dy[ n ];
		}
	}
	return val;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// F   = tensor with Ndim dimensions
// idx = (x,y,..) with Ndim integer elements.
//
// Job of this function is to figure out value of F at
//  index idx. Usually this is just F( idx ).
//
// But when extrapolating,
//  this function needs to figure out value as some kind of linear
//  combination of values where idx is in range.
//
// Strategy: handle each dimension of the idx one by one through
//  a recursion. which_dim tells us where we are in the recursion, from
//  1 to N (and when we hit N+1 we're done).
//
template< typename T, numeric::Size N >
inline
T
get_val(
	MathNTensor< T, N > const & F,
	utility::fixedsizearray1< int, N > idx,
	utility::fixedsizearray1< CatmullRomSplineBoundaryType, N > const & boundary,
	Size const which_dim )
{
	if ( which_dim <= N ) {
		int const Nmax = F.n_bins( which_dim );
		switch ( boundary[ which_dim ] ) {
		case PERIODIC :
			idx[ which_dim ] = utility::modulo( idx[ which_dim ], Nmax ); // wrap
			return get_val( F, idx, boundary, which_dim + 1 );
		case FLAT :
			if ( idx[ which_dim ] < 0 )  idx[ which_dim ] = 0;
			if ( idx[ which_dim ] > Nmax-1 ) idx[ which_dim ] = Nmax-1;
			return get_val( F, idx, boundary, which_dim + 1 );
		case LINEAR :
			if ( idx[ which_dim ] <= 0 ) {
				utility::fixedsizearray1< int, N > idx_0( idx ); idx_0[ which_dim ] = 0;
				utility::fixedsizearray1< int, N > idx_1( idx ); idx_1[ which_dim ] = 1;
				return ( 1 - idx[ which_dim ] ) * get_val( F, idx_0, boundary, which_dim + 1 ) +
					( idx[ which_dim ] )     * get_val( F, idx_1, boundary, which_dim + 1 );
			}
			if ( idx[ which_dim ] > Nmax-1 ) {
				utility::fixedsizearray1< int, N > idx_Nminus1( idx ); idx_Nminus1[ which_dim ] = Nmax - 1;
				utility::fixedsizearray1< int, N > idx_Nminus2( idx ); idx_Nminus2[ which_dim ] = Nmax - 2;
				return ( idx[ which_dim ] - Nmax + 2 ) * get_val( F, idx_Nminus1, boundary, which_dim + 1 ) +
					( Nmax - 1 - idx[ which_dim ] ) * get_val( F, idx_Nminus2, boundary, which_dim + 1 );
			}
			return get_val( F, idx, boundary, which_dim+1 );
			// case CUBIC:
			// This is the specification in Keys, 1981 [from an old MATLAB implementation of
			//  mine... note: was 1-indexed, but should be 0-indexed to match F]
			// Need to instead work out cubic function of idx, like LINEAR above.
			//         if ( idx < 2 )
			//             idx = 1;
			//             F0 = F(3,:) - 3*F(2,:) + 3*F(1,:);
			//             F_patch = [ F0; F(1:3,:)];
			//         elseif ( idx >= N - 1 )
			//             idx = N - 1;
			//             FNplus1 = 3*F(N,:)- 3*F(N-1,:)+F(N-2,:);
			//             F_patch = [ F(N-2:N,:); FNplus1 ];
			//         else
			//             F_patch = F( [idx-1: idx+2], : );
			//         end
		default :
			utility_exit_with_message( "Unrecognized boundary condition!" );
		}
	}
	utility::fixedsizearray1< Size, N > idx_size;
	for ( Size n = 1; n <= N; n++ ) idx_size[n] = Size( idx[n] );
	runtime_assert( F.is_valid_position( idx_size ) );
	Real const & val = F( idx_size );
	return val;
}


// @brief pull out a 4 x 4 ... x 4 patch surrounding desired gridpoint.
// @details
//   handles boundary conditions too, using a kind of recursion.
//
//   some of that extrapolation could be made faster by caching values.
//
//   could be easily generalized to grab 2 N vals around point, e.g.
//    if we want more accurate interpolant
//
// @author rhiju
template< typename T, numeric::Size N >
inline
void
get_patch(
	utility::vector1< T > & F_patch,
	utility::vector1< Real > & d,
	MathNTensor< T, N > const & F,
	utility::fixedsizearray1< Real, N > const & minval,
	utility::fixedsizearray1< Real, N > const & binwidth,
	utility::fixedsizearray1< Real, N > const & xs,
	utility::fixedsizearray1< CatmullRomSplineBoundaryType, N > const & boundary )
{
	utility::fixedsizearray1< int, N > idx;
	for ( Size n = 1; n <= N; n++ ) {
		Real const idx_real = ( xs[n] - minval[n] ) / binwidth[ n ];
		idx[ n ] = std::floor( idx_real );
		d[ n ] = idx_real - idx[ n ];
	}

	// iterate through every point in the 4 x 4 ... 4 F_patch
	Size const four_to_the_n( static_cast<Size>( pow(4, N) ) );
	runtime_assert( F_patch.size() == four_to_the_n );
	for ( Size k = 0; k < four_to_the_n; k++ ) {
		utility::fixedsizearray1< Size, N > idx_F_patch;
		Size c( k );
		// 'decode' which grid point we're at.
		// For example, in 2 dimensions, the last point will be 15 (decimal), which is
		// 33 (base 4) --> (3,3)
		for ( Size n = N; n >= 1; n-- ) {
			idx_F_patch[ n ] = c % 4;
			c /= 4;
		}
		utility::fixedsizearray1< int, N > idx_F;
		for ( Size n = 1; n <= N; n++ ) {
			idx_F[ n ] = idx[ n ] + ( int( idx_F_patch[ n ] ) - 1 ); // convert 0,1,2,3 to -1,0,1,2
		}
		// F_patch is 1-indexed vector1 to match fixedsizearray1.
		// Note however that MathNTensor F is zero-indexed due to historical reasons.
		F_patch[ k+1 ] = get_val( F, idx_F, boundary, 1 );
	}
}

/// @brief Catmull-Rom spline interpolation of an N-dimensional tensor defined on a equispaced grid.
template< typename T, numeric::Size N >
inline
T
polycubic_interpolate_catmull_rom(
	MathNTensor< T, N > const & F,
	utility::fixedsizearray1< Real, N > const & minval,
	utility::fixedsizearray1< Real, N > const & binwidth,
	utility::fixedsizearray1< Real, N > const & xs,
	utility::fixedsizearray1< CatmullRomSplineBoundaryType, N > const & boundary,
	utility::fixedsizearray1< Real, N > & deriv,
	bool const & compute_deriv = true )
{
	// how far are we nudged into central cube, in each dimension 1,2,...N?
	utility::vector1< T > d( N, 0.0 );
	// To define spline, we need a 4 x 4 x... 4 grid of points that surround x.
	utility::vector1< T > F_patch( static_cast<Size>( pow( 4, N ) ), 0.0 );
	get_patch( F_patch, d, F, minval, binwidth, xs, boundary );

	utility::vector1< T >  deriv_vector( N, 0.0 );
	Real const val = catmull_rom_interpolate( F_patch, d, deriv_vector, compute_deriv );
	if ( compute_deriv ) {
		for ( Size n = 1; n <= N; n++ ) deriv[ n ] = deriv_vector[ n ]/ binwidth[ n ];
	}
	return val;
}

/// @brief Catmull-Rom spline interpolation of an N-dimensional tensor defined on a equispaced grid, without derivative computation.
template< typename T, numeric::Size N >
inline
T
polycubic_interpolate_catmull_rom(
	MathNTensor< T, N > const & F,
	utility::fixedsizearray1< Real, N > const & minval,
	utility::fixedsizearray1< Real, N > const & binwidth,
	utility::fixedsizearray1< Real, N > const & xs,
	utility::fixedsizearray1< CatmullRomSplineBoundaryType, N > const & boundary )
{
	utility::fixedsizearray1< Real, N > deriv_dummy;
	return polycubic_interpolate_catmull_rom( F, minval, binwidth, xs, boundary, deriv_dummy, false /* compute_deriv */ );
}


}
}

#endif
