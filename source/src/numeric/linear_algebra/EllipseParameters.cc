// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/linear_algebra/EllipseParameters.cc
/// @brief Container class for ellipse parameters
/// @author Rebecca Alford (rfalford12@gmail.com)

#include <numeric/linear_algebra/EllipseParameters.hh>
#include <numeric/MathMatrix.hh>
#include <numeric/types.hh>

namespace numeric {
namespace linear_algebra {

EllipseParameters::EllipseParameters():
	utility::pointer::ReferenceCount(),
	h_(0),
	k_(0),
	a_(1),
	b_(1),
	rotation_()
{
	numeric::MathMatrix< numeric::Real > rot_matrix( 2, 2 );
	rot_matrix(0,0) = 1;
	rot_matrix(0,1) = 0;
	rot_matrix(1,0) = 0;
	rot_matrix(1,1) = 1;
	rotation_ = rot_matrix;
}

EllipseParameters::EllipseParameters(
	Real h,
	Real k,
	Real a,
	Real b,
	MathMatrix< Real > r
) : utility::pointer::ReferenceCount(),
	h_( h ),
	k_( k ),
	a_( a ),
	b_( b ),
	rotation_( r )
{}

EllipseParameters::~EllipseParameters(){}

EllipseParameters::EllipseParameters( EllipseParameters const & src ) :
	h_( src.h_ ),
	k_( src.k_ ),
	a_( src.a_ ),
	b_( src.b_ ),
	rotation_( src.rotation_ )
{}

EllipseParametersOP
EllipseParameters::clone() const {
	return EllipseParametersOP( new EllipseParameters( *this ) );
}

// add a buffer to ellipse parameters
void EllipseParameters::add_buffer( Real major_buffer, Real minor_buffer ) {
	a_ = a_ + major_buffer;
	b_ = b_ + minor_buffer;
}

// read-only methods
Real EllipseParameters::center_h() const {
	return h_;
}

Real EllipseParameters::center_k() const {
	return k_;
}

Real EllipseParameters::major_radius() const {
	return a_;
}

Real EllipseParameters::minor_radius() const {
	return b_;
}

MathMatrix< Real > EllipseParameters::rotation() const {
	return rotation_;
}

} //numeric
} //linear_algebra






