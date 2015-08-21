// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   kinematic_closure_helpers.hh
/// @brief  Header file for kinematic_closure_helpers.cc
/// @author Daniel J. Mandell

#ifndef INCLUDED_numeric_kinematic_closure_kinematic_closure_helpers_hh
#define INCLUDED_numeric_kinematic_closure_kinematic_closure_helpers_hh

// Rosetta Headers
#include <numeric/types.hh>

// Utility headers

#include <utility/vector1.fwd.hh>


namespace numeric {
namespace kinematic_closure {

// headers
void printVector(const utility::vector1<numeric::Real>& V);
void printMatrix(const utility::vector1<utility::vector1<numeric::Real> >& M);
void multMatrix(const utility::vector1<utility::vector1<numeric::Real> >& A,
	const utility::vector1<utility::vector1<numeric::Real> >& B,
	utility::vector1<utility::vector1<numeric::Real> >& C);
void multTransMatrix(const utility::vector1<utility::vector1<numeric::Real> >& A,
	const utility::vector1<utility::vector1<numeric::Real> >& B,
	utility::vector1<utility::vector1<numeric::Real> >& C);
bool vectorsEqual(const utility::vector1<numeric::Real>& A, const utility::vector1<numeric::Real>& B);

} // end namespace kinematic_closure
} // end namespace numeric

#endif
