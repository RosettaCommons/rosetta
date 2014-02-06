// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/kinematic_closure/vector.cc
/// @author Kale Kundert

// Unit Headers
#include <numeric/types.hh>
#include <numeric/kinematic_closure/vector.hh>
#include <numeric/xyzVector.hh>

// C++ Headers
#include <cmath>
#include <iostream>

namespace numeric {
namespace kinematic_closure {

// Vector Arithmetic Functions

Real dot (Coordinate const &a, Coordinate const &b) { // {{{1
	return a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
}

Coordinate cross (Coordinate const &a, Coordinate const &b) { // {{{1
	Coordinate result (3);

	result[1] = a[2]*b[3] - a[3]*b[2];
	result[2] = a[3]*b[1] - a[1]*b[3];
	result[3] = a[1]*b[2] - a[2]*b[1];

	return result;
}

Coordinate norm (Coordinate const &a) { // {{{1
	Coordinate result (3);

	Real magnitude = sqrt(
			a[1]*a[1] + a[2]*a[2] + a[3]*a[3]);

	result[1] = a[1] / magnitude;
	result[2] = a[2] / magnitude;
	result[3] = a[3] / magnitude;

	return result;
}
// }}}1

// Overloaded Operators
 
std::ostream& operator << (std::ostream &out, ParameterList const &x) { // {{{1
	for (Size i = 1; i <= x.size(); i++) {
		std::string prefix = (i == 1) ? "[" : " ";
		std::string suffix = (i == x.size()) ? "]" : ",";

		out << prefix << x[i] << suffix;
	}
	return out;
}

std::ostream& operator << (std::ostream &out, ParameterMatrix const &xx) { // {{{1
	for (Size i = 1; i <= xx.size(); i++) {
		std::string prefix = (i == 1) ? "[" : " ";
		std::string suffix = (i == xx.size()) ? "]" : ",";

		out << prefix << xx[i] << suffix << std::endl;
	}
	return out;
}

Coordinate operator + (Coordinate const &a, Coordinate const &b) { // {{{1
	Coordinate result (3);

	result[1] = a[1] + b[1];
	result[2] = a[2] + b[2];
	result[3] = a[3] + b[3];

	return result;
}

Coordinate operator - (Coordinate const &a, Coordinate const &b) { // {{{1
	Coordinate result (3);

	result[1] = a[1] - b[1];
	result[2] = a[2] - b[2];
	result[3] = a[3] - b[3];

	return result;
}

Coordinate operator * (Coordinate const &a, Real const &k) { // {{{1
	Coordinate result (3);

	result[1] = a[1] * k;
	result[2] = a[2] * k;
	result[3] = a[3] * k;

	return result;
}

Coordinate operator * (Real const &k, Coordinate const &a) { // {{{1
	Coordinate result (3);

	result[1] = k * a[1];
	result[2] = k * a[2];
	result[3] = k * a[3];

	return result;
}

Coordinate operator / (Coordinate const &a, Real const &k) { // {{{1
	Coordinate result (3);

	result[1] = a[1] / k;
	result[2] = a[2] / k;
	result[3] = a[3] / k;

	return result;
}

Coordinate operator / (Real const &k, Coordinate const &a) { // {{{1
	Coordinate result (3);

	result[1] = k / a[1];
	result[2] = k / a[2];
	result[3] = k / a[3];

	return result;
}
// }}}1

} // end namespace kinematic_closure
} // end namespace numeric
