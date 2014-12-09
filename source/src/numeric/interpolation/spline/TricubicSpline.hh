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
/// Tricubic spline for smoothly interpolating a function in 3 dimensions
///
/// @detailed
///
///
/// @references
/// Numerical Recipes in c++ 2nd edition
/// Ralf Mueller
///
///
/// @authors Steven Combs, Ralf Mueller, Jens Meiler
/// ported to Rosetta by Andrew Leaver-Fay
///
/// @last_modified March 28 2012
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_numeric_interpolation_spline_TricubicSpline_hh
#define INCLUDED_numeric_interpolation_spline_TricubicSpline_hh

#include <numeric/types.hh>
#include <numeric/interpolation/spline/Cubic_spline.fwd.hh>
#include <numeric/MathTensor.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/MathVector.hh>

#include <utility>

namespace numeric {
namespace interpolation {
namespace spline {

class TricubicSpline
{
public:

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	/// construct generic BicubicSpline
	TricubicSpline(){}

	/// copy constructor
	TricubicSpline* Clone() const
	{
	  return new TricubicSpline( *this);
	}

	/////////////////
	// data access //
	/////////////////


	/// get the second order derivatives of the spline
	MathTensor< Real> const & get_dsecox() const
	{
	  return dsecox_;
	}

	MathTensor< Real> const & get_dsecoy() const
	{
	  return dsecoy_;
	}

	MathTensor< Real> const & get_dsecoz() const
	{
	  return dsecoz_;
	}

	MathTensor< Real> const & get_dsecoxy() const
	{
	  return dsecoxy_;
	}

	MathTensor< Real> const & get_dsecoxz() const
	{
	  return dsecoxz_;
	}

	MathTensor< Real> const & get_dsecoyz() const
	{
	  return dsecoyz_;
	}

	MathTensor< Real> const & get_dsecoxyz() const
	{
	  return dsecoxyz_;
	}

	////////////////
	// operations //
	////////////////

	/// @return value at (x, y)
	Real F( Real x, Real y, Real z ) const;

	/// @return partial derivative at (x, y, z) for x
	Real dFdx( Real x, Real y, Real z ) const;

	/// @return partial derivative at (x, y, z) for y
	Real dFdy( Real x, Real y, Real z ) const;

	/// @return partial derivative at (x, y, z) for z
	Real dFdz( Real x, Real y, Real z ) const;

	/// @return value and derivative at (x, y)
	//void FdF( Real x, Real y, Real z, Real & val, Real & dvaldx, Real & dvaldy, Real & dvaldz) const;

	/// train TricubicSpline
	void train
	(
		const BorderFlag BORDER[3],
		const double START[3],
		const double DELTA[3],
		const MathTensor< Real > &RESULTS,
		const bool LINCONT[3],
		const std::pair< Real, Real > FIRSTBE[3]
	);


private:
	BorderFlag border_[3];   ///< controls the behavior at x/y_0 and x/y_dim-1

	Real start_[3];
	Real delta_[3];    ///< gives the arguments as a sequence of equidistant points

	MathTensor< Real> values_;     ///< f(x,y,z)
	MathTensor< Real> dsecox_;     ///< second order derivatives for x -- d**2/dx**2 f(x,y,z)
	MathTensor< Real> dsecoy_;     ///< second order derivatives for y
	MathTensor< Real> dsecoxy_;    ///< second order derivatives for x and y
	MathTensor< Real> dsecoz_;     ///< second order derivatives for z
	MathTensor< Real> dsecoxz_;    ///< second order derivatives for xz
	MathTensor< Real> dsecoyz_;    ///< second order derivatives for yz
	MathTensor< Real> dsecoxyz_;   ///< second order derivatives for x y and z

	std::pair< Real, Real> firstbe_[3]; ///< first order derivative at x_0/dim-1, y_0/dim-1, z_0/dim-1 can be set for BorderFlag FIRSTDER

	bool LinCont_[3];    ///< if the argument x is outside the range decide if the spline should be continued linearly


};

}//end namespace spline
}//end namespace interpolation
}//end namespace numeric




#endif /* BICUBIC_SPLINE_HH_ */
