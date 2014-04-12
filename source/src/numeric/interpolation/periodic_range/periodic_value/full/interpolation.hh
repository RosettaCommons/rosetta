// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/interpolation/periodic_range/periodic_value/full/interpolation.hh
/// @brief  Interpolation of periodic values over periodic range points at full bin width multiples
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @remarks
///  @li For interpolation into arrays with values at { 0, w, 2w, ... } for bin width w
///  @li Bins are numbered from 0,...,n-1 and range from [ (bin)w, (bin+1)w ]
///  @li Bin number the lower index in the interpolated arrays
///
///  array values at     0         w        2w      ...     (n-2)w    (n-1)w
///  bin number          |    0    |    1    |      ...       |   n-1   |
///  array index         0         1         2      ...      n-2       n-1


#ifndef INCLUDED_numeric_interpolation_periodic_range_periodic_value_full_interpolation_hh
#define INCLUDED_numeric_interpolation_periodic_range_periodic_value_full_interpolation_hh


// Package headers
#include <numeric/numeric.functions.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/interpolation/interpolation.hh>

// C++ headers
#include <cassert>
#include <cmath>


namespace numeric {
namespace interpolation {
namespace periodic_range {
namespace periodic_value {
namespace full {


/// @brief Periodic interpolation bin number of a value
template< typename X >
inline
platform::SSize
bin(
	X const & x, // Independent axis value
	X const & w, // Bin width
	platform::SSize const n // Number of bins
)
{
	assert( w > X( 0.0 ) );
	assert( n > 0 );
	platform::SSize const i( static_cast< platform::SSize >( std::floor( x / w ) ) ); // Non-periodic bin number
	return numeric::modulo( i, n );
}


/// @brief Periodic interpolation bin number of a value
template< typename X >
inline
platform::SSize
bin(
	X const & x, // Independent axis value
	X const & w, // Bin width
	platform::SSize const n, // Number of bins
	X & a // Alpha fraction: ( x - x_l ) / ( x_u - x_l ) for bin [ x_l, x_u ]
)
{
	assert( w > X( 0.0 ) );
	assert( n > 0 );
	X const r( x / w );
	platform::SSize const i( static_cast< platform::SSize >( std::floor( r ) ) ); // Non-periodic bin number
	a = r - i;
	assert( ( a >= X( 0.0 ) ) && ( a <= X( 1.0 ) ) );
	return numeric::modulo( i, n );
}


/// @brief Nearest periodic value
template< typename F >
inline
F
nearest_periodic_value(
	F const & a, // Value
	F const & b, // Base value
	F const & p  // Period
)
{
	return a - ( nearest_ssize( ( a - b ) / p ) * p );
}


/// @brief Periodic linearly interpolated value
template< typename X, typename F, template< typename > class A >
inline
F
interpolated(
	X const & x, // Independent axis value
	X const & w, // Bin width
	platform::SSize const n, // Number of bins
	A< F > const & f, // Interpolation array
	F const & p // Function value period
)
{
	assert( n > 0 );
	X a; // Alpha fraction: ( x - x_l ) / ( x_u - x_l ) for bin [ x_l, x_u ]
	platform::SSize const l( bin( x, w, n, a ) );
	assert( ( l >= 0 ) && ( l < n ) );
	platform::SSize const u( numeric::modulo( l + 1, n ) );
	F const & fl( f( l ) );
	return numeric::remainder( numeric::interpolation::interpolated( a, fl, nearest_periodic_value( f( u ), fl, p ) ), p );
}


/// @brief Periodic linearly interpolated value given the bin and alpha fraction
template< typename X, typename F, template< typename > class A >
inline
F
interpolated(
	platform::SSize const l, // Bin number (== lower index)
	X const & a, // Alpha fraction: ( x - x_l ) / ( x_u - x_l ) for bin [ x_l, x_u ]
	platform::SSize const n, // Number of bins
	A< F > const & f, // Interpolation array
	F const & p // Function value period
)
{
	assert( ( l >= 0 ) && ( l < n ) );
	assert( ( a >= X( 0.0 ) ) && ( a <= X( 1.0 ) ) );
	assert( n > 0 );
	platform::SSize const u( numeric::modulo( l + 1, n ) );
	F const & fl( f( l ) );
	return numeric::remainder( numeric::interpolation::interpolated( a, fl, nearest_periodic_value( f( u ), fl, p ) ), p );
}


/// @brief Periodic bilinearly interpolated value
template< typename X, typename F, template< typename > class A >
inline
F
bilinearly_interpolated(
	X const & x1, // Independent axis 1 value
	X const & x2, // Independent axis 2 value
	X const & w, // Bin width
	platform::SSize const n, // Number of bins
	A< F > const & f, // Interpolation array
	F const & p // Function value period
)
{
	assert( w > X( 0.0 ) );
	assert( n > 0 );
	X a1, a2; // Alpha fractions: ( x - x_l ) / ( x_u - x_l ) for bin [ x_l, x_u ]
	platform::SSize const l1( bin( x1, w, n, a1 ) );
	platform::SSize const l2( bin( x2, w, n, a2 ) );
	assert( ( l1 >= 0 ) && ( l1 < n ) );
	assert( ( l2 >= 0 ) && ( l2 < n ) );
	X const b1( X( 1.0 ) - a1 ); // 1 - a1
	X const b2( X( 1.0 ) - a2 ); // 1 - a2
	platform::SSize const u1( numeric::modulo( l1 + 1, n ) );
	platform::SSize const u2( numeric::modulo( l2 + 1, n ) );
	F const fll( f( l1, l2 ) );
	F const ful( nearest_periodic_value( f( u1, l2 ), fll, p ) );
	F const flu( nearest_periodic_value( f( l1, u2 ), fll, p ) );
	F const fuu( nearest_periodic_value( f( u1, u2 ), fll, p ) );
	return numeric::remainder(
	 ( b1 * b2 * fll ) +
	 ( a1 * b2 * ful ) +
	 ( b1 * a2 * flu ) +
	 ( a1 * a2 * fuu ),
	 p );
}


/// @brief Periodic bilinearly interpolated value given the bins and alpha fractions
template< typename X, typename F, template< typename > class A >
inline
F
bilinearly_interpolated(
	platform::SSize const l1, // Axis 1 bin number (== lower index)
	platform::SSize const l2, // Axis 2 bin number (== lower index)
	X const & a1, // Axis 1 alpha fraction: ( x1 - x1_l ) / ( x1_u - x1_l ) for bin [ x1_l, x1_u ]
	X const & a2, // Axis 2 alpha fraction: ( x2 - x2_l ) / ( x2_u - x2_l ) for bin [ x2_l, x2_u ]
	platform::SSize const n, // Number of bins
	A< F > const & f, // Interpolation array
	F const & p // Function value period
)
{
	assert( ( l1 >= 0 ) && ( l1 < n ) );
	assert( ( l2 >= 0 ) && ( l2 < n ) );
	assert( ( a1 >= X( 0.0 ) ) && ( a1 <= X( 1.0 ) ) );
	assert( ( a2 >= X( 0.0 ) ) && ( a2 <= X( 1.0 ) ) );
	assert( n > 0 );
	X const b1( X( 1.0 ) - a1 ); // 1 - a1
	X const b2( X( 1.0 ) - a2 ); // 1 - a2
	platform::SSize const u1( numeric::modulo( l1 + 1, n ) );
	platform::SSize const u2( numeric::modulo( l2 + 1, n ) );
	F const fll( f( l1, l2 ) );
	F const ful( nearest_periodic_value( f( u1, l2 ), fll, p ) );
	F const flu( nearest_periodic_value( f( l1, u2 ), fll, p ) );
	F const fuu( nearest_periodic_value( f( u1, u2 ), fll, p ) );
	return numeric::remainder(
	 ( b1 * b2 * fll ) +
	 ( a1 * b2 * ful ) +
	 ( b1 * a2 * flu ) +
	 ( a1 * a2 * fuu ),
	 p );
}


/// @brief Periodic bilinearly interpolated value and derivatives
template< typename X, typename F, template< typename > class A >
inline
F
bilinearly_interpolated(
	X const & x1, // Independent axis 1 value
	X const & x2, // Independent axis 2 value
	X const & w, // Bin width
	platform::SSize const n, // Number of bins
	A< F > const & f, // Interpolation array
	F const & p, // Function value period
	F & df_dx1, // Derivate wrt axis 1
	F & df_dx2 // Derivate wrt axis 2
)
{
	assert( w > X( 0.0 ) );
	assert( n > 0 );
	X a1, a2; // Alpha fractions: ( x - x_l ) / ( x_u - x_l ) for bin [ x_l, x_u ]
	platform::SSize const l1( bin( x1, w, n, a1 ) );
	platform::SSize const l2( bin( x2, w, n, a2 ) );
	assert( ( l1 >= 0 ) && ( l1 < n ) );
	assert( ( l2 >= 0 ) && ( l2 < n ) );
	X const b1( X( 1.0 ) - a1 ); // 1 - a1
	X const b2( X( 1.0 ) - a2 ); // 1 - a2
	platform::SSize const u1( numeric::modulo( l1 + 1, n ) );
	platform::SSize const u2( numeric::modulo( l2 + 1, n ) );
	F const fll( f( l1, l2 ) );
	F const ful( nearest_periodic_value( f( u1, l2 ), fll, p ) );
	F const flu( nearest_periodic_value( f( l1, u2 ), fll, p ) );
	F const fuu( nearest_periodic_value( f( u1, u2 ), fll, p ) );
	df_dx1 = ( ( b2 * ( ful - fll ) ) + ( a2 * ( fuu - flu ) ) ) / w;
	df_dx2 = ( ( b1 * ( flu - fll ) ) + ( a1 * ( fuu - ful ) ) ) / w;
	return numeric::remainder(
	 ( b1 * b2 * fll ) +
	 ( a1 * b2 * ful ) +
	 ( b1 * a2 * flu ) +
	 ( a1 * a2 * fuu ),
	 p );
}


/// @brief Periodic bilinearly interpolated value and derivatives given the bins and alpha fractions
template< typename X, typename F, template< typename > class A >
inline
F
bilinearly_interpolated(
	platform::SSize const l1, // Axis 1 bin number (== lower index)
	platform::SSize const l2, // Axis 2 bin number (== lower index)
	X const & a1, // Axis 1 alpha fraction: ( x1 - x1_l ) / ( x1_u - x1_l ) for bin [ x1_l, x1_u ]
	X const & a2, // Axis 2 alpha fraction: ( x2 - x2_l ) / ( x2_u - x2_l ) for bin [ x2_l, x2_u ]
	X const & w, // Bin width
	platform::SSize const n, // Number of bins
	A< F > const & f, // Interpolation array
	F const & p, // Function value period
	F & df_dx1, // Derivate wrt axis 1
	F & df_dx2 // Derivate wrt axis 2
)
{
	assert( ( l1 >= 0 ) && ( l1 < n ) );
	assert( ( l2 >= 0 ) && ( l2 < n ) );
	assert( ( a1 >= X( 0.0 ) ) && ( a1 <= X( 1.0 ) ) );
	assert( ( a2 >= X( 0.0 ) ) && ( a2 <= X( 1.0 ) ) );
	assert( w > X( 0.0 ) );
	assert( n > 0 );
	X const b1( X( 1.0 ) - a1 ); // 1 - a1
	X const b2( X( 1.0 ) - a2 ); // 1 - a2
	platform::SSize const u1( numeric::modulo( l1 + 1, n ) );
	platform::SSize const u2( numeric::modulo( l2 + 1, n ) );
	F const fll( f( l1, l2 ) );
	F const ful( nearest_periodic_value( f( u1, l2 ), fll, p ) );
	F const flu( nearest_periodic_value( f( l1, u2 ), fll, p ) );
	F const fuu( nearest_periodic_value( f( u1, u2 ), fll, p ) );
	df_dx1 = ( ( b2 * ( ful - fll ) ) + ( a2 * ( fuu - flu ) ) ) / w;
	df_dx2 = ( ( b1 * ( flu - fll ) ) + ( a1 * ( fuu - ful ) ) ) / w;
	return numeric::remainder(
	 ( b1 * b2 * fll ) +
	 ( a1 * b2 * ful ) +
	 ( b1 * a2 * flu ) +
	 ( a1 * a2 * fuu ),
	 p );
}


} // namespace full
} // namespace periodic_value
} // namespace periodic_range
} // namespace interpolation
} // namespace numeric


#endif // INCLUDED_numeric_interpolation_periodic_range_periodic_value_full_interpolation_HH
