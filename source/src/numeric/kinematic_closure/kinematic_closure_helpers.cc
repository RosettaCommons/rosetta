// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   kinematic_closure_helpers.cc
/// @brief  helpers for functions associated with bridgeObjects.cc
/// @author Daniel J. Mandell

// C++ headers
// Utility headers
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>
// Rosetta Headers
#include <numeric/types.hh>
#include <iostream>

#define SMALL 0.00001  // small epsilon to test for equality

using numeric::Real;

namespace numeric {
namespace kinematic_closure {

// tests if all elements of two vectors are equal
bool vectorsEqual(const utility::vector1<Real>& A, const utility::vector1<Real>& B) {
	for ( unsigned i=1; i<=A.size(); i++ ) {
		if ( std::abs(long(A[i]) - long(B[i])) > SMALL ) {
			return false;
		}
	}
	return true;
}

} // end namespace kinematic_closure
} // end namespace numeric

