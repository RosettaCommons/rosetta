// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


//////////////////////////////////////////////////////////////////////
/// @begin Bicubic_spline
///
/// @brief
/// Bicubic spline for all your hearts desires
///
/// @detailed
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
/// @authors Steven Combs, Ralf Mueller, Jens Meiler
///
/// @last_modified August 20 2010
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_numeric_interpolation_spline_Bicubic_spline_hh
#define INCLUDED_numeric_interpolation_spline_Bicubic_spline_hh

#include <numeric/types.hh>
#include <numeric/interpolation/spline/Cubic_spline.fwd.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/MathVector.hh>

namespace numeric {
namespace interpolation {
namespace spline {

class BicubicSpline
{
public:

   //////////////////////////////////
   // construction and destruction //
	//////////////////////////////////

	/// construct generic BicubicSpline
	BicubicSpline()
	{
	}

	/// copy constructor
	BicubicSpline* Clone() const
	{
		return new BicubicSpline( *this);
	}

   /////////////////
   // data access //
   /////////////////


	/// get the second order derivatives of the spline
	MathMatrix< Real> const &get_dsecox() const
	{
		return dsecox_;
	}

	MathMatrix< Real> const &get_dsecoy() const
	{
		return dsecoy_;
	}

	MathMatrix< Real> const &get_dsecoxy() const
	{
		return dsecoxy_;
	}

	//////////////
	// operator //
	//////////////

	////////////////
	// operations //
	////////////////

	/// @return value at (x, y)
	Real F( const MathVector< Real> &ARGUMENTS) const;

	/// @return partial derivative at (x, y) for x
	Real dFdx( const MathVector< Real> &ARGUMENTS) const;

	/// @return partial derivative at (x, y) for y
	Real dFdy( const MathVector< Real> &ARGUMENTS) const;

	/// @return value at (x, y)
	Real F( Real x, Real y ) const;

	/// @return partial derivative at (x, y) for x
	Real dFdx(  Real x, Real y ) const;

	/// @return partial derivative at (x, y) for y
	Real dFdy(  Real x, Real y ) const;

	/// @return value and derivative at (x, y)
	std::pair< Real, MathVector< Real> > FdF( const MathVector< Real> &ARGUMENTS) const;

	/// train BicubicSpline
	void train (
		const BorderFlag BORDER[2],
		const Real START[2],
		const Real DELTA[2],
		const MathMatrix< Real> &RESULTS,
		const bool LINCONT[2],
		const std::pair< Real, Real> FIRSTBE[2]
	);

	inline  Real sqr ( const Real & x ) const{
		return x*x;
	}


private:
	BorderFlag border_[2];   ///< controls the behavior at x/y_0 and x/y_dim-1

	Real start_[2], delta_[2];    ///< gives the arguments as a sequence of equidistant points

	MathMatrix< Real> values_;     ///< f(x)
	MathMatrix< Real> dsecox_;     ///< second order derivatives for x
	MathMatrix< Real> dsecoy_;     ///< second order derivatives for y
	MathMatrix< Real> dsecoxy_;    ///< second order derivatives for x and y

	std::pair< Real, Real> firstbe_[2]; ///< first order derivative at x_0/dim-1, y_0/dim-1, z_0/dim-1 can be set for BorderFlag FIRSTDER

	bool LinCont_[2];    ///< if the argument x is outside the range decide if the spline should be continued linearly


};

}//end namespace spline
}//end namespace interpolation
}//end namespace numeric




#endif /* BICUBIC_SPLINE_HH_ */
