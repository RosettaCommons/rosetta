// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/angle.functions.hh
/// @brief  Trigonometric functions
/// @author Frank M. D'Ippolito (Objexx@objexx.com)
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_numeric_angle_functions_hh
#define INCLUDED_numeric_angle_functions_hh


// Package headers
#include <numeric/numeric.functions.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/Quaternion.hh>
#include <numeric/xyzMatrix.hh>

// C++ headers
#include <cmath>


namespace numeric {


/// @brief Principal value of angle in radians on ( -pi, pi ]
template< typename T >
inline
T
principal_angle( T const & angle )
{
	return remainder( angle, numeric::NumericTraits< T >::pi_2() );
}


/// @brief Principal value of angle in radians on ( -pi, pi ]
template< typename T >
inline
T
principal_angle_radians( T const & angle )
{
	return remainder( angle, numeric::NumericTraits< T >::pi_2() );
}


/// @brief Principal value of angle in degrees on ( -180, 180 ]
template< typename T >
inline
T
principal_angle_degrees( T const & angle )
{
	return remainder( angle, T( 360.0 ) );
}


/// @brief Positive principal value of angle in radians on [ 0, 2*pi )
template< typename T >
inline
T
nonnegative_principal_angle( T const & angle )
{
	return modulo( angle, numeric::NumericTraits< T >::pi_2() );
}


/// @brief Positive principal value of angle in radians on [ 0, 2*pi )
template< typename T >
inline
T
nonnegative_principal_angle_radians( T const & angle )
{
	return modulo( angle, numeric::NumericTraits< T >::pi_2() );
}


/// @brief Positive principal value of angle in degrees on [ 0, 360 )
template< typename T >
inline
T
nonnegative_principal_angle_degrees( T const & angle )
{
	return modulo( angle, T( 360.0 ) );
}


/// @brief Nearest periodic value of angle to a base angle in radians
template< typename T >
inline
T
nearest_angle( T const & angle, T const & base_angle )
{
	return angle - ( nearest_ssize( ( angle - base_angle ) / numeric::NumericTraits< T >::pi_2() ) * numeric::NumericTraits< T >::pi_2() );
}


/// @brief Nearest periodic value of angle to a base angle in radians
template< typename T >
inline
T
nearest_angle_radians( T const & angle, T const & base_angle )
{
	return angle - ( nearest_ssize( ( angle - base_angle ) / numeric::NumericTraits< T >::pi_2() ) * numeric::NumericTraits< T >::pi_2() );
}


/// @brief Nearest periodic value of angle to a base angle in degrees
template< typename T >
inline
T
nearest_angle_degrees( T const & angle, T const & base_angle )
{
	return angle - ( nearest_ssize( ( angle - base_angle ) / T( 360.0 ) ) * T( 360.0 ) );
}


/// @brief Interconvert Quaternion <=> Rotation Matrix
template <typename T>
inline
void R2quat( xyzMatrix< T > const &R, Quaternion< T > &Q ) {
	T S;
	if ( R.xx() > R.yy() && R.xx() > R.zz() )  {
		S  = sqrt( 1.0 + R.xx() - R.yy() - R.zz() ) * 2;
		Q.x(0.25 * S);
		Q.y((R.xy() + R.yx() ) / S);
		Q.z((R.zx() + R.xz() ) / S);
		Q.w((R.zy() - R.yz() ) / S);
	} else if ( R.yy() > R.zz() ) {
		S  = sqrt( 1.0 + R.yy() - R.xx() - R.zz() ) * 2;
		Q.x((R.yx() + R.xy() ) / S);
		Q.y(0.25 * S);
		Q.z((R.zy() + R.yz() ) / S);
		Q.w((R.xz() - R.zx() ) / S);
	} else {
		S  = sqrt( 1.0 + R.zz() - R.xx() - R.yy() ) * 2;
		Q.x((R.xz() + R.zx() ) / S);
		Q.y((R.zy() + R.yz() ) / S);
		Q.z(0.25 * S);
		Q.w((R.yx() - R.xy()) / S);
	}
}


/// @brief Interconvert Quaternion <=> Rotation Matrix
template <typename T>
inline
void
quat2R( Quaternion< T > const &Q ,xyzMatrix< T > &R ) {
	T xx = Q.x()*Q.x(); T xy = Q.x()*Q.y(); T xz = Q.x()*Q.z();
	T xw = Q.x()*Q.w(); T yy = Q.y()*Q.y(); T yz = Q.y()*Q.z();
	T yw = Q.y()*Q.w(); T zz = Q.z()*Q.z(); T zw = Q.z()*Q.w();

	R.xx(1-2*(yy+zz)); R.xy(  2*(xy-zw)); R.xz(  2*(xz+yw));
	R.yx(  2*(xy+zw)); R.yy(1-2*(xx+zz)); R.yz(  2*(yz-xw));
	R.zx(  2*(xz-yw)); R.zy(  2*(yz+xw)); R.zz(1-2*(xx+yy));
}


} // namespace numeric


#endif // INCLUDED_numeric_angle_functions_HH
