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
/// Base class for abstract N-dimensional PolycubicSpline.
///
/// @author Vikram K. Mulligan (vmullig@uw.edu)
///
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_numeric_interpolation_spline_PolycubicSplineBase_hh
#define INCLUDED_numeric_interpolation_spline_PolycubicSplineBase_hh

#include <numeric/types.hh>
#include <numeric/interpolation/spline/PolycubicSpline.fwd.hh>
#include <numeric/interpolation/spline/PolycubicSplineBase.fwd.hh>
#include <numeric/MathNTensor.hh>

#include <utility>

namespace numeric {
namespace interpolation {
namespace spline {

class PolycubicSplineBase
{
public:

	typedef numeric::Size Size;

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	/// @brief Construct generic PolycubicSplineBase
	///
	PolycubicSplineBase() :
		dimensionality_(0)
	{}

	/// @brief Initialization constructor.
	///
	PolycubicSplineBase( Size const dimensionality_in ) :
		dimensionality_( dimensionality_in )
	{}

	/// @brief Virtual destructor needed for polymorphism.
	///
	virtual ~PolycubicSplineBase() {}

	/////////////////
	// data access //
	/////////////////

	/// @brief Get the dimensionality of the derived class.
	/// @details Needed to allow pointer casts.
	Size dimensionality() const { return dimensionality_; }

protected:

	/// @brief Allow derived class to set the dimensionality.
	///
	void set_dimensionality( Size const dimensionality_in ) { dimensionality_ = dimensionality_in; }

private:

	/// @brief The dimensionality of the derived class.
	///
	Size dimensionality_ = 0;

};


}//end namespace spline
}//end namespace interpolation
}//end namespace numeric


#endif /* INCLUDED_numeric_interpolation_spline_PolycubicSplineBase_hh */
