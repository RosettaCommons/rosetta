// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/interpolation/interpolation.hh
/// @brief  Interpolation functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_interpolation_interpolation_hh
#define INCLUDED_numeric_interpolation_interpolation_hh


// Package headers
#include <numeric/numeric.functions.hh>
#include <numeric/NumericTraits.hh>

// C++ headers
#include <utility/assert.hh>
#include <cmath>


namespace numeric {
namespace interpolation {


/// @brief Linearly interpolated value: f( x )
/// @note  Extrapolates if x not in [ x1, x2 ]
template< typename X, typename F >
inline
F
interpolated(
	X const & x,
	X const & x1,
	X const & x2,
	F const & f1,
	F const & f2
)
{
	assert( x2 - x1 != X( 0 ) );
	return f1 + ( ( x - x1 ) / ( x2 - x1 ) ) * ( f2 - f1 ); // f( x )
}


/// @brief Linearly interpolated value: f( x )
/// @note  Extrapolates if a not in [ 0, 1 ]
template< typename X, typename F >
inline
F
interpolated(
	X const & a, // Alpha fraction: ( x - x1 ) / ( x2 - x1 )
	F const & f1,
	F const & f2
)
{
	return f1 + ( a * ( f2 - f1 ) ); // f( x )
}


/// @brief Linearly interpolated delta value: f( x ) - f1
/// @note  Extrapolates if a not in [ 0, 1 ]
template< typename X, typename F >
inline
F
interpolated_delta(
	X const & a, // Alpha fraction: ( x - x1 ) / ( x2 - x1 )
	F const & f1,
	F const & f2
)
{
	return a * ( f2 - f1 ); // f( x )
}


/// @brief Bilinearly interpolated value: f( x, y )
template< typename X, typename Y, typename F >
inline
F
bilinearly_interpolated(
	X const & x,
	X const & x1,
	X const & x2,
	Y const & y,
	Y const & y1,
	Y const & y2,
	F const & f11, // f( x1, y1 )
	F const & f12, // f( x1, y2 )
	F const & f21, // f( x2, y1 )
	F const & f22  // f( x2, y2 )
)
{
	assert( x2 - x1 != X( 0 ) );
	assert( y2 - y1 != Y( 0 ) );
	X const ax( ( x - x2 ) / ( x2 - x1 ) ); // alpha_x fraction
	Y const ay( ( y - y2 ) / ( y2 - y1 ) ); // alpha_y fraction
	X const bx( X( 1.0 ) - ax ); // beta_x == 1 - alpha_x
	Y const by( Y( 1.0 ) - ay ); // beta_y == 1 - alpha_y
	return
	 ( bx * by * f11 ) +
	 ( bx * ay * f12 ) +
	 ( ax * by * f21 ) +
	 ( ax * ay * f22 );
}


/// @brief Bilinearly interpolated value
template< typename X, typename Y, typename F >
inline
F
bilinearly_interpolated(
	X const & ax,  // alpha_x fraction: ( x - x1 ) / ( x2 - x1 )
	Y const & ay,  // alpha_y fraction: ( y - y1 ) / ( y2 - y1 )
	F const & f11, // f( x1, y1 )
	F const & f12, // f( x1, y2 )
	F const & f21, // f( x2, y1 )
	F const & f22  // f( x2, y2 )
)
{
	X const bx( X( 1.0 ) - ax ); // beta_x == 1 - alpha_x
	Y const by( Y( 1.0 ) - ay ); // beta_y == 1 - alpha_y
	return
	 ( bx * by * f11 ) +
	 ( bx * ay * f12 ) +
	 ( ax * by * f21 ) +
	 ( ax * ay * f22 );
}


/// @brief Bilinearly interpolated value
template< typename X, typename Y, typename F >
inline
F
bilinearly_interpolated(
	X const & ax,  // alpha_x fraction: ( x - x1 ) / ( x2 - x1 )
	Y const & ay,  // alpha_y fraction: ( y - y1 ) / ( y2 - y1 )
	X const & bx,  // beta_x == 1 - alpha_x
	Y const & by,  // beta_y == 1 - alpha_y
	F const & f11, // f( x1, y1 )
	F const & f12, // f( x1, y2 )
	F const & f21, // f( x2, y1 )
	F const & f22  // f( x2, y2 )
)
{
	assert( eq_tol( bx, X( 1.0 ) - ax, NumericTraits< X >::tolerance() * 1000, NumericTraits< X >::tolerance() * 1000 ) );
	assert( eq_tol( by, Y( 1.0 ) - ay, NumericTraits< Y >::tolerance() * 1000, NumericTraits< Y >::tolerance() * 1000 ) );
	return
	 ( bx * by * f11 ) +
	 ( bx * ay * f12 ) +
	 ( ax * by * f21 ) +
	 ( ax * ay * f22 );
}


} // namespace interpolation
} // namespace numeric


#endif // INCLUDED_numeric_interpolation_interpolation_HH
