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

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#endif // SERIALIZATION

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



#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
numeric::SplineParameters::save( Archive & arc ) const {
	arc( CEREAL_NVP( ylo ) ); // platform::Real
	arc( CEREAL_NVP( yhi ) ); // platform::Real
	arc( CEREAL_NVP( y2lo ) ); // platform::Real
	arc( CEREAL_NVP( y2hi ) ); // platform::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
numeric::SplineParameters::load( Archive & arc ) {
	arc( ylo ); // platform::Real
	arc( yhi ); // platform::Real
	arc( y2lo ); // platform::Real
	arc( y2hi ); // platform::Real
}

SAVE_AND_LOAD_SERIALIZABLE( numeric::SplineParameters );

/// @brief Automatically generated serialization method
template< class Archive >
void
numeric::CubicPolynomial::save( Archive & arc ) const {
	arc( CEREAL_NVP( c0 ) ); // platform::Real
	arc( CEREAL_NVP( c1 ) ); // platform::Real
	arc( CEREAL_NVP( c2 ) ); // platform::Real
	arc( CEREAL_NVP( c3 ) ); // platform::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
numeric::CubicPolynomial::load( Archive & arc ) {
	arc( c0 ); // platform::Real
	arc( c1 ); // platform::Real
	arc( c2 ); // platform::Real
	arc( c3 ); // platform::Real
}

SAVE_AND_LOAD_SERIALIZABLE( numeric::CubicPolynomial );
#endif // SERIALIZATION
