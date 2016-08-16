// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   numeric/kinematic_closure/vector.hh
/// @brief  Implement the primitive vector operations (again).
/// @author Kale Kundert
///
/// This file should really not exist, because rosetta already has an xyzVector
/// class which implements all of this functionality.  However, for historical
/// reasons, the kinematic closure algorithms do not use xyzVector and must
/// therefore reimplement the primitive vector arithmetic operations.  Although
/// a better solution would be to refactor the kinematic closure algorithms,
/// this would be a much larger undertaking.

#ifndef INCLUDED_numeric_kinematic_closure_vector_HH
#define INCLUDED_numeric_kinematic_closure_vector_HH

// Unit Headers
#include <numeric/kinematic_closure/types.hh>
#include <numeric/xyzVector.hh>

// C++ Headers
#include <iosfwd>

namespace numeric {
namespace kinematic_closure {

// Vector Arithmetic Functions

Real dot (Coordinate const &a, Coordinate const &b);
Coordinate cross (Coordinate const &a, Coordinate const &b);
Coordinate norm (Coordinate const &a);

// Overloaded Operators

std::ostream& operator << (std::ostream &out, ParameterList const &a);
std::ostream& operator << (std::ostream &out, ParameterMatrix const &a);
Coordinate operator + (Coordinate const &a, Coordinate const &b);
Coordinate operator - (Coordinate const &a, Coordinate const &b);
Coordinate operator * (Coordinate const &a, Real const &k);
Coordinate operator * (Real const &k, Coordinate const &a);
Coordinate operator / (Coordinate const &a, Real const &k);
Coordinate operator / (Real const &k, Coordinate const &a);

// Pseudo-Assignment Operators

// I decided to overload the stream operators to copy utility::vector1 objects
// into xyzVectors and vice versa.  The assignment operator would have been a
// better choice, but it can only be overloaded by member functions and I don't
// want to make any changes to the utility::vector1 class.  Regardless, this is
// really an abuse of the stream operator, so be careful and make sure you know
// what's going on when using it.

template <class T>
Coordinate& operator << (Coordinate &a, xyzVector<T> const &b) {
	a.resize(3);

	a[1] = b.x();
	a[2] = b.y();
	a[3] = b.z();

	return a;
}

template <class T>
xyzVector<T>& operator << (xyzVector<T> &a, Coordinate const &b) {
	a.x() = b[1];
	a.y() = b[2];
	a.z() = b[3];

	return a;
}

} // end namespace kinematic_closure
} // end namespace numeric

#endif

