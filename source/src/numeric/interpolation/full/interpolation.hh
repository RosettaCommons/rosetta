// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/interpolation/full/interpolation.hh
/// @brief  Interpolation over points at full bin width multiples
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @remarks
///  @li Interpolation into arrays with values at i*w for array index i and bin width w


#ifndef INCLUDED_numeric_interpolation_full_interpolation_hh
#define INCLUDED_numeric_interpolation_full_interpolation_hh


// Package headers
#include <numeric/interpolation/interpolation.hh>

// Platform headers
#include <platform/types.hh>

// C++ headers
#include <cmath>


namespace numeric {
namespace interpolation {
namespace full {


namespace bin_width {


/// @brief Lower array index of interpolation bin for an independent axis value
template< typename X >
inline
platform::SSize
l(
	X const & x, // Independent axis value
	X const & w, // Bin width
	X & a // Alpha fraction: ( x - x_l ) / ( x_u - x_l ) for bin [ x_l, x_u ]
)
{
	assert( w > X( 0.0 ) );
	X const r( x / w );
	platform::SSize const l( static_cast< platform::SSize > ( std::floor( r )) );
	a = r - l;
	assert( ( a >= X( 0.0 ) ) && ( a <= X( 1.0 ) ) );
	return l;
}


/// @brief Linearly interpolated value from array
template< typename X, typename F, template< typename > class A >
inline
F
interpolated(
	X const & x, // Independent axis value
	X const & w, // Bin width
	A< F > const & f // Interpolation array
)
{
	X a; // Alpha fraction: ( x - x_l ) / ( x_u - x_l ) for bin [ x_l, x_u ]
	platform::SSize const el( l( x, w, a ) );
	return numeric::interpolation::interpolated( a, f( el ), f( el + 1 ) );
}


} // namespace bin_width


namespace bin_density {


/// @brief Lower array index of interpolation bin for an independent axis value
template< typename X >
inline
platform::SSize
l(
	X const & x, // Independent axis value
	X const & p, // Bins per x unit (inverse bin width)
	X & a // Alpha fraction: ( x - x_l ) / ( x_u - x_l ) for bin [ x_l, x_u ]
)
{
	assert( p > X( 0.0 ) );
	X const r( x * p );
	platform::SSize const l( std::floor( r ) );
	a = r - l;
	assert( ( a >= X( 0.0 ) ) && ( a <= X( 1.0 ) ) );
	return l;
}


/// @brief Linearly interpolated value from array
template< typename X, typename F, template< typename > class A >
inline
F
interpolated(
	X const & x, // Independent axis value
	X const & p, // Bins per x unit (inverse bin width)
	A< F > const & f // Interpolation array
)
{
	X a; // Alpha fraction: ( x - x_l ) / ( x_u - x_l ) for bin [ x_l, x_u ]
	platform::SSize const l( l( x, p, a ) );
	return numeric::interpolation::interpolated( a, f( l ), f( l + 1 ) );
}


} // namespace bin_density


} // namespace full
} // namespace interpolation
} // namespace numeric


#endif // INCLUDED_numeric_interpolation_full_interpolation_HH
