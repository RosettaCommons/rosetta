// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//////////////////////////////////////////////////////////////////////
///
/// @brief
/// read the header file!
///
/// @references
/// Numerical Recipes in c++ 2nd edition
/// Ralf Mueller
///
///
/// @author Steven Combs, Ralf Mueller, Jens Meiler
///
/////////////////////////////////////////////////////////////////////////
#include <numeric/interpolation/spline/Cubic_spline.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/MathMatrix_operations.hh>

#include <numeric/types.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/MathVector.srlz.hh>
#endif // SERIALIZATION


namespace numeric {
namespace interpolation {
namespace spline {


//  train CubicSpline
//  BORDER determines the behavior of the spline at the borders (natural, first derivative, periodic)
//  START the start of the interval the spline is defined on
//  DELTA the distance between two support points of the spline
//  RESULTS the function values at the support points of the spline
//  FIRSTBE values for the first order derivative at begin and end of spline (only FIRSTDER)
CubicSpline &CubicSpline::train
(
	const BorderFlag BORDER,
	const Real START,
	const Real DELTA,
	const MathVector< Real> &RESULTS,
	const std::pair< Real, Real> &FIRSTBE
)
{

	// determine value for dimension of x
	const int dim( RESULTS.size());

	// assigning values
	border_ = BORDER;
	start_  = START;
	delta_  = DELTA;

	// auxiliary variables
	const Real delta1( DELTA / 6), delta2( 2 * DELTA / 3);

	values_  = RESULTS;


	MathMatrix< Real> coeffs( dim, dim);
	MathVector< Real> derivs( dim);

	// train once for the values of fxx
	// those values are equivalent for every type of spline considered here
	for ( int i( 1); i < dim - 1; ++i ) {
		coeffs( i, i - 1) = delta1;
		coeffs( i, i    ) = delta2;
		coeffs( i, i + 1) = delta1;

		derivs( i) = ( values_( i + 1) - 2 * values_( i) + values_( i-1)) / DELTA;
	}

	// setting the second order derivative on start and end to 0, "natural" cubic spline
	if ( border_ == e_Natural ) {
		coeffs(     0,     0) = 1;
		coeffs( dim-1, dim-1) = 1;
	}

	// periodic, the function starts over after reaching the ends, continuously differentiable everywhere
	if ( border_ == e_Periodic ) {
		coeffs(     0, dim - 1) = delta1;
		coeffs(     0,       0) = delta2;
		coeffs(     0,       1) = delta1;
		derivs(     0)          = ( values_( 1) - 2 * values_( 0) + values_( dim-1)) / DELTA;

		coeffs( dim - 1, dim - 2) = delta1;
		coeffs( dim - 1, dim - 1) = delta2;
		coeffs( dim - 1,       0) = delta1;
		derivs( dim - 1)          = ( values_( 0) - 2 * values_( dim-1) + values_( dim-2)) / DELTA;
	}

	// set the first order derivative at x_0 to first_start and at x_dim-1 to first_end
	if ( border_ == e_FirstDer ) {
		coeffs(       0,       0) = -delta2/2;
		coeffs(       0,       1) = -delta1;
		derivs(       0)          = FIRSTBE.first - ( values_( 1) - values_( 0)) / DELTA;

		coeffs( dim - 1, dim - 1) = delta2 / 2;
		coeffs( dim - 1, dim - 2) = delta1;
		derivs( dim - 1)          = FIRSTBE.second - ( values_( dim - 1) - values_( dim - 2)) / DELTA;
	}

	// computation of the second order derivatives in every given point of the spline
	derivs = coeffs.inverse() * derivs;
	dsecox_ = derivs;

	return *this;
}


/// @brief return value at certain ARGUMENT
/// @param ARGUMENT x value
/// @return function value at ARGUMENT
Real CubicSpline::F( const Real &ARGUMENT) const
{
	// number of grid points
	const int dim( values_.size());

	// determine i with start_+(i-1)*delta_ < ARGUMENT < start_+i*delta_ for the correct supporting points
	int i( int( floor( ( ARGUMENT - start_) / delta_)) + 1);

	// not within supporting points - left
	if ( i < 1 ) {
		// if the spline is periodic, adjust i to be positive ( > 0) and within range
		if ( border_ == e_Periodic ) {
			// bring close to actual range
			i %= dim;

			// if between end and start
			if ( i == 0 ) {
				// see Numerical recipes in C++, pages 116-118
				Real dxp( fmod( ARGUMENT - start_, delta_) / delta_); // relative offset from the beginning of the actual interval
				if ( dxp < 0.0 ) {
					dxp += 1.0;
				}

				// return
				return Function( dim - 1, 0, dxp);
			}
			// generate positive index value ( > 0)
			while ( i < 1 )
					{
				i += dim;
			}
		} else {
			return Function( 0, 1, Real( 0)) + ( ARGUMENT - start_) * Derivative( 0, 1, Real( 0));
		}
	} else if ( i >= dim ) {
		// not within supporting points - right
		const int end( dim - 1);

		// if the spline is periodic, adjust i to be positive ( > 0)
		if ( border_ == e_Periodic ) {
			// generate index value within range
			i %= dim;

			// special case, where interpolation happens between end and beginning
			if ( i == 0 ) {
				// see Numerical recipes in C++, pages 116-118
				const Real dxp( fmod( ARGUMENT - start_, delta_) / delta_); // relative offset from the beginning of the actual interval

				// return
				return Function( end, 0, dxp);
			}
		} else {
			// derivative at last supporting points
			return Function( end -1, end, 1.0) + ( ARGUMENT - start_ - ( dim - 1) * delta_) * Derivative( end -1, end, 1.0);
		}
	}

	// see Numerical recipes in C++, pages 116-118
	Real dxp( fmod( ARGUMENT - start_, delta_) / delta_); // delta_akt is an offset from the beginning of the actual interval\n";
	if ( dxp < 0.0 ) {
		dxp += 1.0;
	}

	return Function( i - 1, i, dxp);
}


/// @brief return derivative at certain ARGUMENT
/// @param ARGUMENT x value
/// @return derivative at ARGUMENT
Real CubicSpline::dF( const Real &ARGUMENT) const
{
	// number of grid points
	const int dim( values_.size());

	// determine i with start_+(i-1)*delta_ < ARGUMENT < start_+i*delta_ for the correct supporting points
	int i( int( floor( ( ARGUMENT - start_) / delta_)) + 1);

	// not within supporting points - left
	if ( i < 1 ) {
		// if the spline is periodic, adjust i to be positive ( > 0) and within range
		if ( border_ == e_Periodic ) {
			// bring close to actual range
			i %= dim;

			// if between end and start
			if ( i == 0 ) {
				// see Numerical recipes in C++, pages 116-118
				Real dxp( fmod( ARGUMENT - start_, delta_) / delta_); // relative offset from the beginning of the actual interval
				if ( dxp < 0.0 ) {
					dxp += 1.0;
				}

				// return
				return Derivative( dim - 1, 0, dxp);
			}
			// generate positive index value ( > 0)
			while ( i < 1 )
					{
				i += dim;
			}
		} else {
			return Derivative( 0, 1, Real( 0));
		}
	} else if ( i >= dim ) {
		// not within supporting points - right
		const int end( dim - 1);

		// if the spline is periodic, adjust i to be positive ( > 0)
		if ( border_ == e_Periodic ) {
			// generate index value within range
			i %= dim;

			// special case, where interpolation happens between end and beginning
			if ( i == 0 ) {
				// see Numerical recipes in C++, pages 116-118
				const Real dxp( fmod( ARGUMENT - start_, delta_) / delta_); // relative offset from the beginning of the actual interval

				// return
				return Derivative( end, 0, dxp);
			}
		} else {
			// derivative at last supporting points
			return Derivative( end -1, end, 1.0);
		}
	}

	// see Numerical recipes in C++, pages 116-118
	Real dxp( fmod( ARGUMENT - start_, delta_) / delta_); // delta_akt is an offset from the beginning of the actual interval\n";
	if ( dxp < 0.0 ) {
		dxp += 1.0;
	}

	return Derivative( i - 1, i, dxp);
}

/// @brief return derivative and value at certain ARGUMENT
/// @param ARGUMENT x value
/// @return value and derivative at ARGUMENT
std::pair< Real, Real> CubicSpline::FdF( const Real &ARGUMENT) const
{
	return std::pair<Real, Real>( F( ARGUMENT), dF( ARGUMENT));
}


///////////////////////
// helper  functions //
///////////////////////

/// @brief calculate function between two cells
/// @param INDEX_LEFT index of left grid point
/// @param INDEX_RIGHT index of right grid point
/// @param DXP relative distance from left grid point, must be element [0, 1]
/// @return function depending on relative distance DXP
Real CubicSpline::Function( const int INDEX_LEFT, const int INDEX_RIGHT, const Real DXP) const
{
	// see Numerical recipes in C++, pages 116-118
	// relative distance from right grid point
	const Real dxm( 1 - DXP);
	const Real dx3p( ( DXP * DXP * DXP - DXP) * sqr( delta_) / 6); // =0 at the gridpoints, adds cubic part of the spline
	const Real dx3m( ( dxm * dxm * dxm - dxm) * sqr( delta_) / 6); // =0 at the gridpoints, adds cubic part of the spline

	return
		dxm * values_( INDEX_LEFT) + DXP  * values_( INDEX_RIGHT)
		+ dx3m * dsecox_( INDEX_LEFT) + dx3p * dsecox_( INDEX_RIGHT);
}


/// @brief calculate derivative between two cells
/// @param INDEX_LEFT index of left grid point
/// @param INDEX_RIGHT index of right grid point
/// @param DXP relative distance from left grid point, must be element [0, 1]
/// @return derivative depending on relative distance DXP
Real CubicSpline::Derivative( const int INDEX_LEFT, const int INDEX_RIGHT, const Real DXP) const
{
	// see Numerical recipes in C++, pages 116-118
	// relative distance from right grid point
	const Real dxm( 1 - DXP);

	return
		( values_( INDEX_RIGHT) - values_( INDEX_LEFT)) / delta_
		- ( 3 * dxm * dxm - 1) / 6 * delta_ * dsecox_( INDEX_LEFT)
		+ ( 3 * DXP * DXP - 1) / 6 * delta_ * dsecox_( INDEX_RIGHT);
}

}//spline
}//interpolation
}//numeric


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
numeric::interpolation::spline::CubicSpline::save( Archive & arc ) const {
	arc( CEREAL_NVP( border_ ) ); // enum numeric::interpolation::spline::BorderFlag
	arc( CEREAL_NVP( start_ ) ); // Real
	arc( CEREAL_NVP( delta_ ) ); // Real
	arc( CEREAL_NVP( values_ ) ); // MathVector<Real>
	arc( CEREAL_NVP( dsecox_ ) ); // MathVector<Real>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
numeric::interpolation::spline::CubicSpline::load( Archive & arc ) {
	arc( border_ ); // enum numeric::interpolation::spline::BorderFlag
	arc( start_ ); // Real
	arc( delta_ ); // Real
	arc( values_ ); // MathVector<Real>
	arc( dsecox_ ); // MathVector<Real>
}

SAVE_AND_LOAD_SERIALIZABLE( numeric::interpolation::spline::CubicSpline );
#endif // SERIALIZATION
