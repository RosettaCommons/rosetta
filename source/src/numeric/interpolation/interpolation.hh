// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <utility/fixedsizearray1.hh>
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

/// @brief Perform cubic interpolation over each of N axes, using the
/// 2^N derivatives at 2^N gridpoints
/// @details
/// The way encoding gridpoints and derivatives into a linear structure like this
/// is actually pretty simple. Imagine the "right or left" part of a cube, or the
/// "derivative taken or not" on a particular variable, as zero or one. Then "just
/// the actual function value" maps to 000, the z derivative (for example) maps
/// to 001, d2/dydz maps to 011, etc.
/// @param[in] n_derivs is a 2^N x 2^N matrix: 2^N derivatives at 2^N gridpoints
/// @param[in] dbbp is how far along the bin our target point is, in each direction
/// @param[in] binwbb is the bin width in each direction
/// @param[out] val is the interpolated value
/// @param[out] dvaldbb are the interpolated derivatives
template < Size N >
void
polycubic_interpolation(
	utility::fixedsizearray1< utility::fixedsizearray1< Real, ( 1 << N ) >, ( 1 << N ) > n_derivs,
	utility::fixedsizearray1< Real, N > dbbp,
	utility::fixedsizearray1< Real, N > binwbb,
	Real & val,
	utility::fixedsizearray1< Real, N > & dvaldbb
) {
	utility::fixedsizearray1< Real, N > invbinwbb;
	utility::fixedsizearray1< Real, N > binwbb_over_6;
	utility::fixedsizearray1< Real, N > dbbm;
	utility::fixedsizearray1< Real, N > dbb3p;
	utility::fixedsizearray1< Real, N > dbb3m;
	for ( Size ii = 1; ii <= N; ++ii ) {
		invbinwbb[ ii ] = 1/binwbb[ ii ];
		binwbb_over_6[ ii ] = binwbb[ ii ] / 6 ;
		dbbm[ ii ] = 1 - dbbp[ ii ];
		dbb3p[ ii ] = ( dbbp[ ii ] * dbbp[ ii ] * dbbp[ ii ] - dbbp[ ii ] ) * binwbb[ ii ] * binwbb_over_6[ ii ];
		dbb3m[ ii ] = ( dbbm[ ii ] * dbbm[ ii ] * dbbm[ ii ] - dbbm[ ii ] ) * binwbb[ ii ] * binwbb_over_6[ ii ];
	}

	// there are 2^N deriv terms, i.e. the value, the N first derivatives,
	// the N^2 second derivatives... up to the single Nth derivative

	// The value has its own functional form.
	val = 0;
	for ( Size iid = 1; iid <= (1 << N); ++iid ) {
		for ( Size iiv = 1; iiv <= (1 << N); ++iiv ) {
			Real valterm = n_derivs[ iid ][ iiv ];

			for ( Size jj = 1; jj <= N; ++jj ) { // each bb
				Size two_to_the_jj_compl = 1 << ( N - jj );
				if ( ( iiv - 1 ) & two_to_the_jj_compl ) {
					valterm *= ( ( iid - 1 ) & two_to_the_jj_compl ) ? dbb3p[ jj ] : dbbp[ jj ];
				} else {
					valterm *= ( ( iid - 1 ) & two_to_the_jj_compl ) ? dbb3m[ jj ] : dbbm[ jj ];
				}
			}

			val += valterm;
		}
	}

	//Each of the N first derivatives.
	for ( Size bbn = 1; bbn <= N; ++bbn ) {
		dvaldbb[ bbn ] = 0;

		for ( Size iid = 1; iid <= (1 << N); ++iid ) {
			for ( Size iiv = 1; iiv <= (1 << N); ++iiv ) {
				Real valterm = n_derivs[ iid ][ iiv ]; // v000

				for ( Size jj = 1; jj <= N; ++jj ) {
					Size two_to_the_jj_compl = 1 << ( N - jj );

					// Half of the values from iiv = 1 to 2^N come from
					// "bb_bin_next" and half from "bb_bin."
					if ( ( iiv - 1 ) & two_to_the_jj_compl ) {

						// Does the iid-th derivative have a jj-backbone deriv?
						if ( ( iid - 1 ) & two_to_the_jj_compl ) {
							valterm *= ( bbn == jj ) ?      ( 3 * dbbp[ jj ] * dbbp[ jj ] - 1 ) * binwbb_over_6[ jj ] : dbb3p[ jj ];
						} else {
							valterm *= ( bbn == jj ) ?      invbinwbb[ jj ]                                           :  dbbp[ jj ];
						}
					} else { // bb_bin
						if ( ( iid - 1 ) & two_to_the_jj_compl ) { // is derived
							valterm *= ( bbn == jj ) ? -1 * ( 3 * dbbm[ jj ] * dbbm[ jj ] - 1 ) * binwbb_over_6[ jj ] : dbb3m[ jj ];
						} else {
							// subtract all terms where bbn was taken from bb_bin
							valterm *= ( jj == bbn ) ? -1 * invbinwbb[ jj ]                                           :  dbbm[ jj ];
						}
					}
				}

				dvaldbb[ bbn ] += valterm;
			}
		}
	}
}


} // namespace interpolation
} // namespace numeric


#endif // INCLUDED_numeric_interpolation_interpolation_HH
