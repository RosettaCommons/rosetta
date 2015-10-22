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
#include <basic/basic.hh>
#include <ObjexxFCL/Fmath.hh>    // for mod
#include <cassert>               // for assert
#include <cmath>                 // for fabs, sqrt

// Numeric headers
#include <numeric/constants.hh>  // for pi_2


//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end


// C++ Headers

namespace basic {


//     util_basic.cc - general utility functions that don't fit into any
//     other util*.cc files


////////////////////////////////////////////////////////////////////////////////
///
/// @brief calculates quadratic polynomial solutions
///
/// @details
///
///     solves for solutions of x in the polynomial: a*x^2+b*x+c=0
///
/// @param[in]   a - in - x^2 term
/// @param[in]   b - in - x term
/// @param[in]   c - in - constant term
/// @param[out]   n1 - out - one solution
/// @param[out]   n2 - out - another solution
///
/// @remarks courtesy of Jerry Tsai
///
/// @references
///
/// @author ctsa 8=19-03
///
/////////////////////////////////////////////////////////////////////////////////
void
calc_quadratic(
	double a,
	double b,
	double c,
	double & n1,
	double & n2
)
{
	//cj

	double bsq = b*b;
	double ac4 = 4*a*c;
	double st = std::sqrt( bsq - ac4 );

	//cj    std::cout << F( 8, 3, a ) << ' ' << F( 8, 3, b ) << ' ' << F( 8, 3, c ) << std::endl;
	//cj    std::cout << F( 8, 3, bsq ) << ' ' << F( 8, 3, ac4 ) << ' ' << F( 8, 3, st ) << std::endl;

	n1 = ((-b)+st)/(2*a);
	n2 = ((-b)-st)/(2*a);
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief subtract angles in degrees, restricting the range of the result
///
/// @details
///
///     given angles a and b in degrees, get a-b restricted to
///     [-180.,180.), assuming that a-b=a-b+n*360, n=any integer
///
/// @param[in]   a - in - angle in degrees
/// @param[in]   b - in - angle in degrees
///
/// @return  angle a-b in degrees, restricted to the specified range
///
/// @remarks
///
/// @references
///
/// @author ctsa 8-19-03
///
/////////////////////////////////////////////////////////////////////////////////
double
subtract_degree_angles(
	double a,
	double b
)
{
	return periodic_range( a - b, double(360.0) );
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief subtract angles in radians, restricting the range of the result
///
/// @details
///
///     given angles a and b in degrees, get a-b restricted to
///     [-pi,pi), assuming that a-b=a-b+n*2*pi, n=any integer
///
/// @param[in]   a - in - angle in radians
/// @param[in]   b - in - angle in radians
///
/// @return  angle a-b in degrees, restricted to the specified range
///
/// @remarks
///
/// @references
///
/// @author ctsa 8-19-03
///
/////////////////////////////////////////////////////////////////////////////////
double
subtract_radian_angles(
	double a,
	double b
)
{
	using namespace numeric::constants::f;
	return periodic_range( a - b, pi_2 );
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief  a is restricted to [-x/2,x/2), assuming that a=a+n*x, n=any integer
///
/// @details
///
/// @param[in]   a - in - input value with periodicity x
/// @param[in]   x - in - periodicity of a
///
/// @return  a restricted to [-x/2,x/2)
///
/// @remarks
///
/// @references
///
/// @author ctsa 8-19-03
///
/////////////////////////////////////////////////////////////////////////////////
double
periodic_range(
	double a,
	double x
)
{
	double const halfx = 0.5f * x;
	return ( ( a >= halfx || a < -halfx ) ? mod( mod( a, x ) + ( x + halfx ), x ) - halfx : a );
}

////////////////////////////////////////////////////////////////////////////////
///
///
/// @details
///
/// @param[in]   a - in - input value with periodicity x
/// @param[in]   x - in - periodicity of a
///
/// @return  a restricted to [0.,x)
///
/// @remarks
///
/// @references
///
/// @author ctsa 8-19-03
///
/////////////////////////////////////////////////////////////////////////////////
double
unsigned_periodic_range(
	double a,
	double x
)
{
	return ( ( a >= x || a < 0.0 ) ? mod( mod( a, x ) + x, x ) : a );
}

/// @brief taken from wobble.cc
void
angle_in_range( double & ang )
{
	int const odd = (int)std::fabs( double( mod( static_cast< int >( ang / double(180.0) ), 2 ) ) ); // 1 if ang/180 is odd,0 if ang/180 is even
	ang = mod( ang, double(180.0) );
	if ( odd == 0 ) return;
	if ( ang > 0.0 ) {
		ang -= double(180.0);
		return;
	}
	ang += double(180.0);
	assert( ang <= double(180.0) && ang > double(-180.0));
}


} // namespace basic
