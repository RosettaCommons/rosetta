// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#ifndef INCLUDED_basic_basic_hh
#define INCLUDED_basic_basic_hh


// util_basic Function Declarations

namespace basic {

void
calc_quadratic(
	double a,
	double b,
	double c,
	double & n1,
	double & n2
);


double
subtract_degree_angles(
	double a,
	double b
);


double
subtract_radian_angles(
	double a,
	double b
);


double
periodic_range(
	double a,
	double x
);

/// @brief a is restricted to [0.,x), assuming that a=a+n*x,, n=any integer
double
unsigned_periodic_range(
	double a,
	double x
);


/// taken from wobble.cc
void
angle_in_range( double & ang );


}

#endif
