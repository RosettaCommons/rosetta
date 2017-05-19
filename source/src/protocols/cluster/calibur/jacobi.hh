// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/kalngyk/jacobi.hh
/// @author YK Ng & SC Li (kalngyk@gmail.com)

#ifndef external_calibur_jacobi_HH
#define external_calibur_jacobi_HH

#include <cmath>
#include <iostream>

namespace protocols {
namespace cluster {
namespace calibur {

/**
* Computes the eigenvalues and (normalized) eigenvectors of a real symmetric
* matrix a[3][3] in d[3] and v[3][3] respectively.
*
* Only the upper-right triangle of a[3][3] is referenced and over-written.
*
* (The method follows closely the discussions in "Numerical Recipes in C".)
*/
bool jacobi3( double a[3][3], double d[3], double v[3][3] );

#ifdef _TEST_JACOBI_
void print_matrix(double a[3][3]);
#endif

}
}
}

#endif
