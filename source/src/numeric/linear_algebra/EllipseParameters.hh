// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/linear_algebra/EllipseParameters.hh
/// @brief Container class for ellipse parameters
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_numeric_linear_algebra_EllipseParameters_hh
#define INCLUDED_numeric_linear_algebra_EllipseParameters_hh

#include <numeric/linear_algebra/EllipseParameters.fwd.hh>

// Project Headers
#include <numeric/MathMatrix.hh>
#include <numeric/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace numeric {
namespace linear_algebra {

///@brief Container class for ellipse parameters
class EllipseParameters : public utility::pointer::ReferenceCount {

public:

	/// @brief Default constructor for ellipse parameters
	EllipseParameters();

	/// @Brief Construct a set of ellipse parameters given axes, center, and rotation
	EllipseParameters(
		Real h,
		Real k,
		Real a,
		Real b,
		MathMatrix< Real > rotation
	);

	/// @brief Create a new set of ellipse parameters from an existing initialized object
	EllipseParameters( EllipseParameters const & src );

	virtual ~EllipseParameters();

	EllipseParametersOP
	clone() const;

	/// @brief Add buffer to major & minor radius
	void add_buffer( Real major_buffer, Real minor_buffer );

	/// @brief Get the x-coordinate of the ellipse center
	Real center_h() const;

	/// @brief Get the y-coordinate of the ellipse center
	Real center_k() const;

	/// @brief Get the major radius length of the ellipse
	Real major_radius() const;

	/// @brief Get the minor radius length of the ellipse
	Real minor_radius() const;

	/// @brief 2D matrix defining the ellipse rotation
	MathMatrix< Real > rotation() const;

private:

	// ellipse center (x coord)
	Real h_;

	// ellipse center (y coord)
	Real k_;

	// ellipse major axis
	Real a_;

	// ellipse minor axis
	Real b_;

	// ellipse rotation
	MathMatrix< Real > rotation_;

};


} //numeric
} //linear_algebra



#endif //INCLUDED_numeric_linear_algebra_EllipseParameters_hh





