// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/deriv/distance_deriv.hh
/// @brief  inline function for computing f1/f2 derivatives for a function of an angle
/// @author Phil Bradley did all the hard work deriving the math represented here.
/// @author Andrew Leaver-Fay copy-and-pasted Phil's code into this file from
/// the AngleConstraint.cc file for general use.

#ifndef INCLUDED_numeric_deriv_angle_deriv_hh
#define INCLUDED_numeric_deriv_angle_deriv_hh

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <utility/assert.hh>

namespace numeric {
namespace deriv   {

/////////////////////////////////////////////////////////////////////////////
// accumulate F1, F2 contributions from terms that look like
//
//       d(M-F)
// dot( -------- , w )
//       d phi
//
// where phi is moving M and F is fixed.
//
// F1 collects terms f1 that look like: (-u_phi) dot f1
// F2 collects terms f2 that look like: (-u_phi x R_phi) dot f2
//
// where u_phi is the unit vector axis of phi and R_phi is a point
// on the axis
//
// basic id: d(M-F)/dphi = u_phi x ( M - R_phi )
//
// and dot( a, b x c ) = dot( b, c x a )
//
//
template < class P >
inline
void
helper(
	xyzVector< P > const & M,
	xyzVector< P > const & w,
	xyzVector< P > & F1,
	xyzVector< P > & F2
)
{
	F2 += w;
	F1 += cross( w, M );
}

/////////////////////////////////////////////////////////////////////////////
// calculates f1,f2 contributions for dtheta_dphi
//
// where phi is a torsion angle moving p1 while p2 and p3 are fixed
template < class P >
inline
void
p1_theta_deriv(
	xyzVector< P > const & p1,
	xyzVector< P > const & p2,
	xyzVector< P > const & p3,
	xyzVector< P > & F1,
	xyzVector< P > & F2
)
{
	typedef P Real;
	typedef xyzVector< P > Vector;

	using numeric::conversions::radians;

	// to avoid problems with dtheta/dx around 0 and 180 degrees
	// truncate x a bit in the calculation of the derivative
	static Real const small_angle( radians( Real(0.1) ) );
	static Real const big_angle( radians( Real(179.9) ) );
	static Real const max_x( std::cos( small_angle ));
	static Real const min_x( std::cos( big_angle ));
	// dtheta_dx has a value of ~ 572.96 for min_x and max_x
	// this goes to infinity as x goes to -1 or 1

	Vector v1( p1 - p2 );
	Vector v2( p3 - p2 );
	Real const n1( v1.length() );
	Real const n2( v2.length() );
	if ( n1 < Real(1e-9) || n2 < Real(1e-9) ) {
		return;
	}

	// calculate dx/dphi where x = cos theta = dot( v1,v2) / (n1*n2)

	Real x = dot(v1,v2)/(n1*n2);

	// only v1 depends on phi, not v2 (we are assuming only p1 is moving)
	// so we get two terms, one from v1 on top and one from n1 = |v1| on
	// the bottom

	Vector f1(0.0),f2(0.0);

	{ // first term
		Real const f = Real(1.0) / ( n1 * n2 );
		helper( p1, f * v2, f1, f2 );
	}

	{ // second term
		Real const f = Real(-1.0) * x / ( n1 * n1 );
		helper( p1, f * v1, f1, f2 );
	}

	x = std::min( std::max( min_x, x ), max_x );
	Real const dtheta_dx = Real(-1.0) / sqrt( Real(1.0) - x*x );
	f1 *= dtheta_dx;
	f2 *= dtheta_dx;

	// translation of p1 along v1 or perpendicular to v1 and v2 ==> deriv=0
	assert( f1.distance( cross(f2,p1) ) < Real(1e-3) && // see helper fcn
		std::abs( dot( f2, v1 ) ) < Real(1e-3) &&
		std::abs( dot( f2, cross( v1, v2 ) ) ) < Real(1e-3) );


	{ // more debugging
		// pretend axis = u2, R_phi = p2
		ASSERT_ONLY(Vector const u_phi( v2.normalized() );)
		ASSERT_ONLY(Vector const R_phi( p2 );)
		ASSERT_ONLY(Real const deriv = - dot( u_phi, f1 ) - dot( cross( u_phi, R_phi ), f2);)
		assert( std::abs( deriv ) < Real(1e-3) );
		//std::cout << "deriv: " << deriv<< ' ' <<
		// F(9,3,u_phi(1)) << F(9,3,u_phi(2)) << F(9,3,u_phi(3)) << ' ' <<
		// F(9,3,R_phi(1)) << F(9,3,R_phi(2)) << F(9,3,R_phi(3)) << "\nF1,F2: " <<
		// F(9,3,f1(1)) << F(9,3,f1(2)) << F(9,3,f1(3)) << ' ' <<
		// F(9,3,f2(1)) << F(9,3,f2(2)) << F(9,3,f2(3)) << std::endl;
	}


	F1 += f1;
	F2 += f2;
}

/////////////////////////////////////////////////////////////////////////////
// calculates x and dtheta_dx where
// v1 = p1-p2,
// v2 = p3-p2
// n1 = norm( v1 )
// n2 = norm( v2 )
// x  = dot(v1,v2) / (n1*n2), and
// dtheta_dx = Real(-1.0) / sqrt( Real(1.0) - x*x );
template < class P >
inline
void
x_and_dtheta_dx(
	xyzVector< P > const & v1,
	xyzVector< P > const & v2,
	P const n1,
	P const n2,
	P & x,
	P & dtheta_dx
)
{
	typedef P Real;

	using numeric::conversions::radians;
	static Real const small_angle( radians( Real(0.1) ) );
	static Real const big_angle( radians( Real(179.9) ) );
	static Real const max_x( std::cos( small_angle ));
	static Real const min_x( std::cos( big_angle ));

	x = dot(v1,v2)/(n1*n2);
	x = std::min( std::max( min_x, x ), max_x );
	dtheta_dx = Real(-1.0) / sqrt( Real(1.0) - x*x );

}

/////////////////////////////////////////////////////////////////////////////
// calculates f1,f2 contributions for dtheta_dphi where
// v1 = p1-p2,
// v2 = p3-p2
// n1 = norm( v1 )
// n2 = norm( v2 )
// x  = dot(v1,v2) / (n1*n2), and
// dtheta_dx = Real(-1.0) / sqrt( Real(1.0) - x*x );
//
// where phi is a torsion angle moving p1 while p2 and p3 are fixed
template < class P >
inline
void
p1_theta_deriv(
	xyzVector< P > const & p1,
	xyzVector< P > const & v1,
	xyzVector< P > const & v2,
	P const n1,
	P const n2,
	P const x,
	P const dtheta_dx,
	xyzVector< P > & F1,
	xyzVector< P > & F2
)
{
	typedef P Real;
	typedef xyzVector< P > Vector;

	using numeric::conversions::radians;

	// to avoid problems with dtheta/dx around 0 and 180 degrees
	// truncate x a bit in the calculation of the derivative
	//static Real const small_angle( radians( Real(0.1) ) );
	//static Real const big_angle( radians( Real(179.9) ) );
	// dtheta_dx has a value of ~ 572.96 for min_x and max_x
	// this goes to infinity as x goes to -1 or 1

	//Vector v1( p1 - p2 );
	//Vector v2( p3 - p2 );
	//Real const n1( v1.length() );
	//Real const n2( v2.length() );
	if ( n1 < Real(1e-9) || n2 < Real(1e-9) ) {
		return;
	}

	// calculate dx/dphi where x = cos theta = dot( v1,v2) / (n1*n2)

	//Real x = dot(v1,v2)/(n1*n2);

	// only v1 depends on phi, not v2 (we are assuming only p1 is moving)
	// so we get two terms, one from v1 on top and one from n1 = |v1| on
	// the bottom

	Vector f1(0.0),f2(0.0);

	{ // first term
		Real const f = Real(1.0) / ( n1 * n2 );
		helper( p1, f * v2, f1, f2 );
	}

	{ // second term
		Real const f = Real(-1.0) * x / ( n1 * n1 );
		helper( p1, f * v1, f1, f2 );
	}

	f1 *= dtheta_dx;
	f2 *= dtheta_dx;

	// translation of p1 along v1 or perpendicular to v1 and v2 ==> deriv=0
	//assert( f1.distance( cross(f2,p1) ) < Real(1e-3) && // see helper fcn
	//    std::abs( dot( f2, v1 ) ) < Real(1e-3) &&
	//    std::abs( dot( f2, cross( v1, v2 ) ) ) < Real(1e-3) );


	//{ // more debugging
	// // pretend axis = u2, R_phi = p2
	// Vector const u_phi( v2.normalized() ), R_phi( p2 );
	// Real const deriv = - dot( u_phi, f1 ) - dot( cross( u_phi, R_phi ), f2);
	// assert( std::abs( deriv ) < Real(1e-3) );
	// //std::cout << "deriv: " << deriv<< ' ' <<
	// // F(9,3,u_phi(1)) << F(9,3,u_phi(2)) << F(9,3,u_phi(3)) << ' ' <<
	// // F(9,3,R_phi(1)) << F(9,3,R_phi(2)) << F(9,3,R_phi(3)) << "\nF1,F2: " <<
	// // F(9,3,f1(1)) << F(9,3,f1(2)) << F(9,3,f1(3)) << ' ' <<
	// // F(9,3,f2(1)) << F(9,3,f2(2)) << F(9,3,f2(3)) << std::endl;
	//}


	F1 += f1;
	F2 += f2;
}


/////////////////////////////////////////////////////////////////////////////
/// @brief compute f1/f2 atom derivative vectors for the first point defining
/// an angle. Templated on the precision of the coordinates being represented.
/// Returns the angle, theta, defined by p1->p2->p3, which should be used
/// to evaluate dE_dtheta.  dE_dtheta should then be multiplied into both
/// derivative vectors that are returned. The values of the output variables
/// f1 and f2 are overwritten. Theta is computed in radians.
template < class P >
inline
void
angle_p1_deriv(
	xyzVector< P > const & p1,
	xyzVector< P > const & p2,
	xyzVector< P > const & p3,
	P & theta,
	xyzVector< P > & f1,
	xyzVector< P > & f2
)
{
	//std::cout << "core.mm.MMBondAngleEnergy:     p1 angle deriv:" <<
	//" p1 (" << p1.x() << "," << p1.y() << "," << p1.z() <<
	//") p2 (" << p2.x() << "," << p2.y() << "," << p2.z() <<
	//") p3 (" << p3.x() << "," << p3.y() << "," << p3.z() << ")" <<
	//std::endl;


	typedef P Real;
	typedef xyzVector< P > Vector;

	f1 = Real(0.0);
	f2 = Real(0.0);
	theta = 0.0;

	Vector u1( p1 - p2 );
	Vector u2( p3 - p2 );
	Real const n1_n2( u1.length() * u2.length() );
	if ( n1_n2 < Real(1e-12) ) {
		std::cout << "AngleConstraint::p1_deriv: short bonds: " << n1_n2 <<
			std::endl;
		return;
	}

	p1_theta_deriv( p1, p2, p3, f1, f2 );

	Real d( dot(u1,u2) / n1_n2 );
	Real const tol(1.0e-8);
	if ( d <= Real(-1.0) + tol ) {
		//std::cout << "out-of-tol: " << d << ' ' << std::endl;
		d = Real(-1.0) + tol;
	} else if ( d >= Real(1.0) - tol ) {
		//std::cout << "out-of-tol: " << d << std::endl;
		d = Real(1.0) - tol;
	}
	theta = numeric::arccos( d );

}


/// @brief compute f1/f2 atom derivative vectors for the middle point defining
/// an angle. Templated on the precision of the coordinates being represented.
/// Returns the angle, theta, defined by p1->p2->p3, which should be used
/// to evaluate the dfunc.  dE_dtheta should then be multiplied into both
/// derivative vectors that are returned. The values of the output variables
/// f1 and f2 are overwritten.  Theta is computed in radians.
template < class P >
inline
void
angle_p2_deriv(
	xyzVector< P > const & p1,
	xyzVector< P > const & p2,
	xyzVector< P > const & p3,
	P & theta,
	xyzVector< P > & f1,
	xyzVector< P > & f2
)
{
	//std::cout << "core.mm.MMBondAngleEnergy:     p2 angle deriv:" <<
	//" p1 (" << p1.x() << "," << p1.y() << "," << p1.z() <<
	//") p2 (" << p2.x() << "," << p2.y() << "," << p2.z() <<
	//") p3 (" << p3.x() << "," << p3.y() << "," << p3.z() << ")" <<
	//std::endl;

	typedef P Real;
	typedef xyzVector< P > Vector;

	f1 = Real(0.0);
	f2 = Real(0.0);
	theta = 0.0;

	Vector v1( p1 - p2 );
	Vector v2( p3 - p2 );
	Real const v12( v1.length() * v2.length() );
	if ( v12 < Real(1e-12) ) return;

	// here we use the trick that theta = pi - alpha - beta
	//
	// where alpha and beta are the other angles in the p1,p2,p3 triangle
	// so dtheta_dphi  = -dalpha_dphi - dbeta_dphi
	//
	// for these we can use the p1_theta_deriv formula
	//
	f1 = Real(0.0);
	f2 = Real(0.0);
	p1_theta_deriv( p2, p1, p3, f1, f2 ); // alpha deriv
	p1_theta_deriv( p2, p3, p1, f1, f2 ); // beta deriv

	// translation of p2 atom perpendicular to plane ==> deriv = 0
	//std::cout << "p2 deriv check: " << std::endl;
	assert( std::abs( dot( f2, cross( v1,v2) ) ) < Real(1e-3) );

	Real d( dot(v1,v2) / v12 );
	Real const tol(0.001);
	if ( d <= Real(-1.0) + tol ) {
		//std::cout << "out-of-tol: " << d << ' ' << std::endl;
		d = Real(-1.0) +tol;
	} else if ( d >= Real(1.0) - tol ) {
		//std::cout << "out-of-tol: " << d << std::endl;
		d = Real(1.0) -tol;
	}
	theta = numeric::arccos( d );

	f1 *= Real(-1.0);
	f2 *= Real(-1.0);

}

/// @details useful for computing all three f1/f2 derivative vectors simultaneously.
/// Basically a small refactoring / copy and pasting of the above code.
template < class P >
inline
void
angle_p1_p2_p3_deriv(
	xyzVector< P > const & p1,
	xyzVector< P > const & p2,
	xyzVector< P > const & p3,
	P & theta,
	xyzVector< P > & f1_p1,
	xyzVector< P > & f2_p1,
	xyzVector< P > & f1_p2,
	xyzVector< P > & f2_p2,
	xyzVector< P > & f1_p3,
	xyzVector< P > & f2_p3
)
{
	typedef P Real;
	typedef xyzVector< P > Vector;

	f1_p1 = f2_p1 = f1_p2 = f2_p2 = f1_p3 = f2_p3 = Real(0.0);
	theta = 0.0;

	Vector v12( p1 - p2 );
	Vector v32( p3 - p2 );

	Real const n12( v12.length() );
	Real const n23( v32.length() );

	Real const n12n23( n12 * n23 );
	if ( n12n23 < Real(1e-12) ) return;

	Real x13, dtheta_dx13;
	x_and_dtheta_dx( v12, v32, n12, n23, x13, dtheta_dx13 );
	p1_theta_deriv( p1, v12, v32, n12, n23, x13, dtheta_dx13, f1_p1, f2_p1 );
	p1_theta_deriv( p3, v32, v12, n23, n12, x13, dtheta_dx13, f1_p3, f2_p3 );

	/// Point 2 Derivatives. Comments stolen from above:
	// here we use the trick that theta = pi - alpha - beta
	//
	// where alpha and beta are the other angles in the p1,p2,p3 triangle
	// so dtheta_dphi  = -dalpha_dphi - dbeta_dphi
	//
	// for these we can use the p1_theta_deriv formula


	// alpha deriv for point 2
	Vector v21 = -1 * v12;
	Vector v31 = p3 - p1;
	Real const n13 = v31.length();
	Real alpha_x, alpha_dtheta_dx;
	x_and_dtheta_dx( v21, v31, n12, n13, alpha_x, alpha_dtheta_dx );
	p1_theta_deriv(  p2, v21, v31, n12, n13, alpha_x, alpha_dtheta_dx, f1_p2, f2_p2 );

	// beta deriv for point 2
	Vector v23 = -1 * v32;
	Vector v13 = -1 * v31;
	Real beta_x, beta_dtheta_dx;
	x_and_dtheta_dx( v23, v13, n23, n13, beta_x, beta_dtheta_dx );
	p1_theta_deriv(  p2, v23, v13, n23, n13, beta_x, beta_dtheta_dx, f1_p2, f2_p2 );

	f1_p2 *= Real( -1.0 );
	f2_p2 *= Real( -1.0 );

	Real d( dot(v12,v32) / n12n23 );
	Real const tol(0.001);
	if ( d <= Real(-1.0) + tol ) {
		//std::cout << "out-of-tol: " << d << ' ' << std::endl;
		d = Real(-1.0) +tol;
	} else if ( d >= Real(1.0) - tol ) {
		//std::cout << "out-of-tol: " << d << std::endl;
		d = Real(1.0) -tol;
	}
	theta = numeric::arccos( d );

	/*{ // debug point 1
	Vector test_p1_f1, test_p1_f2; Real test_theta;
	angle_p1_deriv( p1, p2, p3, test_theta, test_p1_f1, test_p1_f2 );
	if ( f1_p1.distance_squared( test_p1_f1 ) > 1e-4 ) {
	std::cout << "Point 1 f1 vector in error: " << f1_p1.x() << " " << f1_p1.y() << " " << f1_p1.z() <<
	" vs " << test_p1_f1.x() << " " << test_p1_f1.y() << " " << test_p1_f1.z() << std::endl;
	}

	if ( f2_p1.distance_squared( test_p1_f2 ) > 1e-4 ) {
	std::cout << "Point 1 f2 vector in error: " << f2_p1.x() << " " << f2_p1.y() << " " << f2_p1.z() <<
	" vs " << test_p1_f2.x() << " " << test_p1_f2.y() << " " << test_p1_f2.z() << std::endl;
	}
	}

	{ // debug point 2
	Vector test_p2_f1, test_p2_f2; Real test_theta;
	angle_p2_deriv( p1, p2, p3, test_theta, test_p2_f1, test_p2_f2 );
	if ( f1_p2.distance_squared( test_p2_f1 ) > 1e-4 ) {
	std::cout << "Point 1 f1 vector in error: " << f1_p2.x() << " " << f1_p2.y() << " " << f1_p2.z() <<
	" vs " << test_p2_f1.x() << " " << test_p2_f1.y() << " " << test_p2_f1.z() << std::endl;
	}

	if ( f2_p2.distance_squared( test_p2_f2 ) > 1e-4 ) {
	std::cout << "Point 1 f2 vector in error: " << f2_p2.x() << " " << f2_p2.y() << " " << f2_p2.z() <<
	" vs " << test_p2_f2.x() << " " << test_p2_f2.y() << " " << test_p2_f2.z() << std::endl;
	}
	}

	{ // debug point 3
	Vector test_p3_f1, test_p3_f2; Real test_theta;
	angle_p1_deriv( p3, p2, p1, test_theta, test_p3_f1, test_p3_f2 );
	if ( f1_p3.distance_squared( test_p3_f1 ) > 1e-4 ) {
	std::cout << "Point 1 f1 vector in error: " << f1_p3.x() << " " << f1_p3.y() << " " << f1_p3.z() <<
	" vs " << test_p3_f1.x() << " " << test_p3_f1.y() << " " << test_p3_f1.z() << std::endl;
	}

	if ( f2_p3.distance_squared( test_p3_f2 ) > 1e-4 ) {
	std::cout << "Point 1 f2 vector in error: " << f2_p3.x() << " " << f2_p3.y() << " " << f2_p3.z() <<
	" vs " << test_p3_f2.x() << " " << test_p3_f2.y() << " " << test_p3_f2.z() << std::endl;
	}
	}*/

}

}
}

#endif


