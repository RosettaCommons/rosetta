// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/cubic_polynomial.hh
/// @brief Functions to evaluate cubic polynomials and cubic polynomial derivatives
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_numeric_cubic_polynomial_hh
#define INCLUDED_numeric_cubic_polynomial_hh

#include <platform/types.hh>

namespace numeric {

struct CubicPolynomial {
	platform::Real c0, c1, c2, c3;
	CubicPolynomial() : c0(0), c1(0), c2(0), c3(0) {}
};


/// @brief %SplineParameters is a simple struct for holding the cubic spline polynomials used in
/// the etable to interpolate the lennard-jones attractive and LK-solvation terms to zero smoothly.
/// These splines have exactly two knots to represent them, and the same x values are used for all
/// the knots: thus the only parameters needed are the y values at the knots, and the second-derivatives
/// for the polynomials at knots.
struct SplineParameters {
	platform::Real ylo,  yhi;
	platform::Real y2lo, y2hi;
	SplineParameters() :
		ylo(0.0),
		yhi(0.0),
		y2lo(0.0),
		y2hi(0.0)
	{}
};

/// @brief Compute cubic polynomial coefficients from a set of SplineParameters
CubicPolynomial
cubic_polynomial_from_spline(
	platform::Real xlo,
	platform::Real xhi,
	SplineParameters const & sp
);


/// @brief Evaluate cubic polynomial at value x given polynomial coefficients
platform::Real
eval_cubic_polynomial(
	platform::Real const x,
	CubicPolynomial const & cp
);

/// @brief Evaluate derivative of cubic polynomial given x and polynomial coefficients
platform::Real
cubic_polynomial_deriv(
	platform::Real const x,
	CubicPolynomial const & cp
);

} // numeric

#endif // numeric_cubic_polynomial_hh

