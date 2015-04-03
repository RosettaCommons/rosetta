// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/deriv/dihedral_deriv.hh
/// @brief  inline function for computing f1/f2 derivatives for a function of a dihedral
/// @author Phil Bradley did all the hard work deriving the math represented here.
/// @author Andrew Leaver-Fay copy-and-pasted Phil's code into this file from
/// the DihedralConstraint.cc file for general use.

#ifndef INCLUDED_numeric_deriv_dihedral_deriv_hh
#define INCLUDED_numeric_deriv_dihedral_deriv_hh

#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/trig.functions.hh>
#include <numeric/conversions.hh>
#include <utility/assert.hh>

#include <float.h>

namespace numeric {
namespace deriv   {

/////////////////////////////////////////////////////////////////////////////
// contributions for a term that looks like
//
//   d(M-F)
// ( ------ x v ) dot w
//   d phi
//
// F1 collects terms f1 that look like: (-u_phi) dot f1
// F2 collects terms f2 that look like: (-u_phi x R_phi) dot f2
//
//
// where u_phi is the torsion axis of rotation
// and R_phi is a terminal atom of this axis
template < class P >
inline
void
helper(
	xyzVector< P > const & M,
	xyzVector< P > const & v,
	xyzVector< P > const & w,
	xyzVector< P > & F1,
	xyzVector< P > & F2
)
{
	typedef xyzVector< P > Vector;

	Vector const f2 = cross(v,w);
	F2 += f2;
	F1 += cross(f2,M);
}

/// @brief Whether computing the f1/f2 derivatives for a function of a dihedral
/// for the end points (p1 or p4) or the middle points (p2 or p3), the second
/// half of the computation is the same.  This "second" function expects
/// the F1 and F2 arrays to have been partially computed, as well as the cosine
/// of the angle theta.  It scales the F1 and F2 vectors by
/// dtheta_dthetaU * dthetaU_dx
template < class P >
inline
void
dihedral_deriv_second(
	xyzVector< P > const & p1,
	xyzVector< P > const & p2,
	xyzVector< P > const & p3,
	xyzVector< P > const & p4,
	P x,
	P & theta,
	xyzVector< P > & F1,
	xyzVector< P > & F2
)
{
	typedef P Real;

	// to avoid problems with dtheta/dx around 0 and 180 degrees
	// truncate x a bit in the calculation of the derivative
	Real const small_angle( numeric::conversions::radians( Real(0.1) ) );
	Real const big_angle( numeric::conversions::radians( Real(179.9) ) );
	Real const max_x( std::cos( small_angle ));
	Real const min_x( std::cos( big_angle ));
	// dtheta_dx has a value of ~ 572.96 for min_x and max_x
	// this goes to infinity as x goes to -1 or 1

	// unsigned version of theta
	ASSERT_ONLY(Real const thetaU( numeric::arccos( x ));)

	theta = dihedral_radians( p1, p2, p3, p4 );

	assert( std::abs( std::abs( theta ) - thetaU ) < 1e-2 );

	x = std::min( std::max( min_x, x ), max_x );
	Real const dthetaU_dx = -1 / sqrt( 1- x*x );
	Real const dtheta_dthetaU( theta < 0 ? -1 : 1 );

	F1 *= dtheta_dthetaU * dthetaU_dx;
	F2 *= dtheta_dthetaU * dthetaU_dx;

}

/// @brief The first half of the computation of the f1/f2 derivatives for
/// an end point of a dihedral.  The values in the output-parameter
/// vectors F1 and F2 are overwritten.  The cosine of the dihedral theta
/// is returned in the output parameter x.
template < class P >
inline
void
dihedral_p1_cosine_deriv_first(
	xyzVector< P > const & p1,
	xyzVector< P > const & p2,
	xyzVector< P > const & p3,
	xyzVector< P > const & p4,
	P & x,
	bool & colinearity_flag,
	xyzVector< P > & F1,
	xyzVector< P > & F2
)
{
	typedef P Real;
	typedef xyzVector< P > Vector;

	F1 = 0.0;
	F2 = 0.0;
	colinearity_flag = false;

	Vector v1( p1-p2 );
	Vector v2( p2-p3 );
	Vector v3( p3-p4 );

	Vector v12( cross( v1, v2 ));
	Vector v23( cross( v2, v3 ));

	Real const n12( v12.length() );
	Real const n23( v23.length() );

	if ( n12 < Real(1e-9) || n23 < Real(1e-9) ) {
		colinearity_flag = true;
		return;
	}

	x = dot( v12, v23 ) / ( n12 * n23 );

	// first term:
	{
		Real const f( Real(1.0) / ( n12 * n23 ) );
		helper( p1, f * v2, v23 , F1, F2);
	}

	// second term
	{
		Real const f( Real(-1.0) * x / ( n12 * n12 ) );
		helper( p1, f * v2, v12, F1, F2 );
	}


	// debugging
	// translation of p1 in the place spanned by v1 and v2
	// does not affect the torsion angle
	// ==> rotation of p1 about an axis perpendicular to this plane
	// also does not change the torsion angle, ie deriv should be 0
	assert( std::abs( dot( F2, v1 ) ) < Real(1e-3) );
	assert( std::abs( dot( F2, v2 ) ) < Real(1e-3) );
	assert( std::abs( dot( F1, cross( v1, v2 ) ) ) < Real(1e-3) );
}

/// @brief compute f1/f2 atom derivative vectors for one of the two end points defining
/// a dihedral angle for some function F. Templated on the precision of the coordinates
/// being represented.  Returns the dihedral angle, theta, defined by p1->p2->p3->p4,
/// which should be used to evaluate the derivative of the function F.  dF_dtheta should
/// then be multiplied into both derivative vectors that are returned. The values of the
/// output variables f1 and f2 are overwritten.  Theta is computed in radians.
/// @brief compute f1/f2 atom derivative vectors for one of the two middle points defining
/// a dihedral angle for some function F. Templated on the precision of the coordinates
/// being represented.  Returns the dihedral angle, theta, defined by p1->p2->p3->p4,
/// which should be used to evaluate the derivative of the function F.  dF_dtheta should
/// then be multiplied into both derivative vectors that are returned. The values of the
/// output variables f1 and f2 are overwritten.  Theta is computed in radians.
template < class P >
inline
void
dihedral_p1_cosine_deriv(
	xyzVector< P > const & p1,
	xyzVector< P > const & p2,
	xyzVector< P > const & p3,
	xyzVector< P > const & p4,
	P & theta,
	xyzVector< P > & f1,
	xyzVector< P > & f2
)
{
	typedef P Real;

	Real x( 0.0 ); bool colinearity;
	dihedral_p1_cosine_deriv_first( p1, p2, p3, p4, x, colinearity, f1, f2 );
	if ( colinearity ) return;
	dihedral_deriv_second( p1, p2, p3, p4, x, theta, f1, f2 );

}

/// @brief The first half of the computation of the f1/f2 derivatives for
/// a central point of a dihedral.  The values in the output-parameter
/// vectors F1 and F2 are overwritten.  The cosine of the dihedral theta
/// is returned in the output parameter x.
template < class P >
inline
void
dihedral_p2_cosine_deriv_first(
	xyzVector< P > const & p1,
	xyzVector< P > const & p2,
	xyzVector< P > const & p3,
	xyzVector< P > const & p4,
	P & x,
	bool & colinearity_flag,
	xyzVector< P > & F1,
	xyzVector< P > & F2
)
{
	typedef P Real;
	typedef xyzVector< P > Vector;

	//std::cout << "p2_cosine_deriv!" << std::endl;

	F1 = Real(0.0);
	F2 = Real(0.0);
	colinearity_flag = false;

	Vector v1( p1-p2 );
	Vector v2( p2-p3 );
	Vector v3( p3-p4 );

	Vector v12( cross( v1, v2 ));
	Vector v23( cross( v2, v3 ));

	Real const n12( v12.length() );
	Real const n23( v23.length() );

	if ( n12 < Real(1e-9) || n23 < Real(1e-9) ) {
		colinearity_flag = true;
		return;
	}

	x = dot( v12, v23) / ( n12 * n23 );

	// here we are taking derivatives of an expression that looks like
	//
	//                   dot( v12, v23 )
	// x = cos theta =  -----------------
	//                      n12 * n23
	//
	// where theta is our dihedral angle


	{ // derivatives of the numerator
		// v1 and v2 both depend on position of p2

		{ // first term
			Real const f( Real(-1.0)/ ( n12 * n23 ) );
			helper( p2, f * v2, v23 , F1, F2);
		}

		{ // second term
			Real const f( Real(-1.0)/ ( n12 * n23 ) );
			helper( p2, f * v1, v23 , F1, F2);
		}

		{ // third term
			Real const f( Real(1.0)/ ( n12 * n23 ) );
			helper( p2, f * v3, v12 , F1, F2);
		}
	}

	{ // derivatives of the denominator
		// v1 and v2 both depend on position of p2

		{ // first term
			Real const f( x / ( n12 * n12 ) );
			helper( p2, f * v2, v12, F1, F2 );
		}

		{ // second term
			Real const f( x / ( n12 * n12 ) );
			helper( p2, f * v1, v12, F1, F2 );
		}

		{ // third term
			Real const f( Real(-1.0) * x / ( n23 * n23 ) );
			helper( p2, f * v3, v23, F1, F2 );
		}
	}

	// debugging
	// translation of p2 along v2 does not change the torsion angle
	assert( std::abs( dot( F2, v2 ) ) < Real(1e-3) );

}


/// @brief compute f1/f2 atom derivative vectors for one of the two middle points defining
/// a dihedral angle for some function F. Templated on the precision of the coordinates
/// being represented.  Returns the dihedral angle, theta, defined by p1->p2->p3->p4,
/// which should be used to evaluate the derivative of the function F.  dF_dtheta should
/// then be multiplied into both derivative vectors that are returned. The values of the
/// output variables f1 and f2 are overwritten.  Theta is computed in radians.
template < class P >
inline
void
dihedral_p2_cosine_deriv(
	xyzVector< P > const & p1,
	xyzVector< P > const & p2,
	xyzVector< P > const & p3,
	xyzVector< P > const & p4,
	P & theta,
	xyzVector< P > & f1,
	xyzVector< P > & f2
)
{
	typedef P Real;

	Real x( 0.0 ); bool colinearity;
	dihedral_p2_cosine_deriv_first( p1, p2, p3, p4, x, colinearity, f1, f2 );
	if ( colinearity ) return;
	dihedral_deriv_second( p1, p2, p3, p4, x, theta, f1, f2 );

}


}
}

#endif


