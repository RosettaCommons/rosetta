// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/cubic_polynomial.cc
/// @brief Functions to evaluate cubic polynomials and cubic polynomial derivatives
/// @author Rebecca Alford (rfalford12@gmail.com)

#include <numeric/cubic_polynomial.hh>

#include <platform/types.hh>

namespace numeric {

CubicPolynomial
cubic_polynomial_from_spline(
	platform::Real xlo,
	platform::Real xhi,
	SplineParameters const & sp
) {
	platform::Real a( xlo ), b( xhi ), c( sp.yhi ), d( sp.ylo ), e( sp.y2hi ), f( sp.y2lo );
	CubicPolynomial cp;
	cp.c0 = ( (b*b*b*f - a*a*a*e)/(b-a) + (a*e - b*f) * (b-a) )  / 6 + (  b*d - a*c ) / ( b-a );
	cp.c1 = ( 3*a*a*e/(b-a) - e*(b-a) + f*(b-a) - 3*b*b*f / (b-a) ) / 6 + ( c - d ) / (b-a);
	cp.c2 = ( 3*b*f - 3*a*e ) / ( 6 * (b-a) );
	cp.c3 = ( e-f ) / ( 6 * (b-a) );
	return cp;
}

platform::Real
eval_cubic_polynomial(
	platform::Real const x,
	CubicPolynomial const & cp
) {
	return ((cp.c3*x+cp.c2)*x+cp.c1)*x + cp.c0;
}

platform::Real
cubic_polynomial_deriv(
	platform::Real const x,
	CubicPolynomial const & cp
) {
	return (3*cp.c3*x + 2*cp.c2)*x + cp.c1;
}

} // numeric


