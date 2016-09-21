// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


//////////////////////////////////////////////////////////////////////
///
/// @brief
/// Cubic spline for all your evil desires
///
/// @details
/// The below comments are for the Bicubic spline but apply for the cubic spline.
/// This is an implementation of an algorithm from Numerical Recipes. It relies heavily
/// on the implementation of the cubic spline (from Numerical Recipes), the MathMatrix, and the
/// MathVector. You MUST USE the MathVector and MathMatrix implementations to use this function.
/// The spline is very customizable and allows you to define the border behaviors (start/end of spline values).
/// If you use the e_Natural (enum) BorderFlag, the start/end (border) of the spline will be linear. This may
/// not be ideal for scoring functions as you want the a smoothing effect at the start/end (border) values.
/// Instead, you probably want to use the e_FirstDeriv (enum) BorderFlag with the first derivate (private member value
/// firstbe_) set to 0. This will cause a smoothing out of the start/end (border) of spline. If you want the splie to
/// be continuous, you should use the e_Periodic (enum) BorderFlag.
///
/// To "train" the spline, use must provide the spline with a MathMatrix (numeric::MathMatrix). Lets look at an example.
///             x values
///    y
///           _1__2__ 3_
///       .1 | 1  2   3 |
///    v  .3 | 4  5   6 |
///    a  .5 | 7  8   9 |
///    l     |__________|
///    u
///    e
///    s
///
/// Given the above Matrix (MathMatrix) You would want your start (START[2] private member value start_) values to be START[] = {1,.1}.
/// You would then want to assign the delta (DELTA[2], private member value delta_) values to DELTA[] = {1,.2}. These delta values
/// is the change between your x values and your y values. For example, the change between x1 and x2 is 1. Therefore, the delta for the
/// x-values will be 1. For y values, you have y.1, y.3 which is a change of .2, therefore the delta will be .2. You do not have to
/// specify an end because the algorithm will stop when it reaches the last value in the matrix.
///
/// Finally, the LinCont determins that if the argument x or y is outside the range decide if the spline should be continued linearly.
///
///
/// @references
/// Numerical Recipes in c++ 2nd edition
/// Ralf Mueller
///
///
/// @author Steven Combs, Ralf Mueller, Jens Meiler
///
/////////////////////////////////////////////////////////////////////////


#ifndef INCLUDED_numeric_interpolation_spline_Cubic_spline_hh
#define INCLUDED_numeric_interpolation_spline_Cubic_spline_hh

#include <numeric/types.hh>
#include <numeric/MathVector.hh>


#include <numeric/interpolation/spline/Cubic_spline.fwd.hh>

namespace numeric {
namespace interpolation {
namespace spline {

class CubicSpline
{
public:

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	//! @brief construct generic CubicSpline
	CubicSpline():
		border_( e_Natural ),
		start_(0.0),
		delta_(0.0)
	{
	}

	//! @brief copy constructor
	CubicSpline* clone() const
	{
		return new CubicSpline( *this);
	}


	CubicSpline & train
	(
		BorderFlag const BORDER,
		Real const START,
		Real const DELTA,
		MathVector< Real> const &RESULTS,
		std::pair< Real, Real> const &FIRSTBE
	);


	Real F( const Real &ARGUMENT) const;

	inline  Real sqr ( const Real x ) const{
		return x*x;
	}


	//! @brief return derivative at ARGUMENT
	//! @param ARGUMENT x value
	//! @return derivative at ARGUMENT
	Real dF( const Real &ARGUMENT) const;

	//! @brief return value and derivative at ARGUMENT
	//! @param ARGUMENT x value
	//! @return value and derivative at ARGUMENT
	std::pair< Real, Real> FdF( const double &ARGUMENT) const;

	//////////////////
	////data access/////
	////////////////////


	//! @brief get the second order derivatives of the spline
	//! @return the second order derivatives at the support points of the spline
	MathVector< Real> const & get_dsecox() const
	{
		return dsecox_;
	}

	//! @brief access to the start value
	//! @return the start of the interval the spline is defined on
	Real get_start() const
	{
		return start_;
	}

	//! @brief access to the delta value
	//! @return the distance between two support points of the spline
	Real get_delta() const
	{
		return delta_;
	}

	//! @brief access to the values
	//! @return the function values at the support points of the spline
	const MathVector< Real> & get_values() const
	{
		return values_;
	}

	bool
	operator == ( CubicSpline const & rhs ) const {
		return border_ == rhs.border_ && start_ == rhs.start_ && delta_ == rhs.delta_ &&
			values_ == rhs.values_ && dsecox_ == rhs.dsecox_;
	}

	bool
	operator != ( CubicSpline const & rhs ) const {
		return ! ( *this == rhs );
	}

private:
	BorderFlag border_; //!< controls the behavior at x_0 and x_dim-1
	Real start_, delta_; //!< gives the arguments as a sequence of equidistant points
	MathVector<Real> values_; //!< f(x)
	MathVector<Real> dsecox_; //!< second order derivatives


	//! @brief calculate function between two cells
	//! @param INDEX_LEFT index of left grid point
	//! @param INDEX_RIGHT index of right grid point
	//! @param DXP relative distance from left grid point, must be element [0, 1]
	//! @return function depending on relative distance DXP
	Real Function( const int INDEX_LEFT, const int INDEX_RIGHT, const Real DXP) const;


	//! @brief calculate derivative between two cells
	//! @param INDEX_LEFT index of left grid point
	//! @param INDEX_RIGHT index of right grid point
	//! @param DXP relative distance from left grid point, must be element [0, 1]
	//! @return derivative depending on relative distance DXP
	Real Derivative( const int INDEX_LEFT, const int INDEX_RIGHT, const Real DXP) const;


#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}//namespace
}
}


#endif
