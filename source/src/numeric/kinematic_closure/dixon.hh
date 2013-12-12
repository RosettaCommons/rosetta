// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   dixon.hh
/// @brief  Header file for dixon code for ALC
/// @author Evangelos A. Coutsias
/// @author Daniel J. Mandell

#ifndef INCLUDED_numeric_kinematic_closure_dixon_hh
#define INCLUDED_numeric_kinematic_closure_dixon_hh

// Rosetta Headers
#include <numeric/types.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>

#include <utility/vector1.fwd.hh>


namespace numeric {
namespace kinematic_closure {

/*
void dixon(
		const utility::vector1<utility::vector1<numeric::Real> >& A,
		const utility::vector1<utility::vector1<numeric::Real> >& B,
		const utility::vector1<utility::vector1<numeric::Real> >& C,
		const utility::vector1<utility::vector1<numeric::Real> >& D,
		const utility::vector1<int>& order,
		utility::vector1<utility::vector1<numeric::Real> >& cos,
		utility::vector1<utility::vector1<numeric::Real> >& sin,
		utility::vector1<utility::vector1<numeric::Real> >& tau,
		int & nsol);
*/

void dixon_eig(
		const utility::vector1<utility::vector1<numeric::Real> >& A,
		const utility::vector1<utility::vector1<numeric::Real> >& B,
		const utility::vector1<utility::vector1<numeric::Real> >& C,
		const utility::vector1<utility::vector1<numeric::Real> >& D,
		const utility::vector1<int>& order,
		utility::vector1<utility::vector1<numeric::Real> >& cos,
		utility::vector1<utility::vector1<numeric::Real> >& sin,
		utility::vector1<utility::vector1<numeric::Real> >& tau,
		int & nsol);

void dixon_sturm(
		const utility::vector1<utility::vector1<numeric::Real> >& A,
		const utility::vector1<utility::vector1<numeric::Real> >& B,
		const utility::vector1<utility::vector1<numeric::Real> >& C,
		const utility::vector1<utility::vector1<numeric::Real> >& D,
		const utility::vector1<int>& order,
		utility::vector1<utility::vector1<numeric::Real> >& cos,
		utility::vector1<utility::vector1<numeric::Real> >& sin,
		utility::vector1<utility::vector1<numeric::Real> >& tau,
		int & nsol);


} // end namespace kinematic_closure
} // end namespace numeric

#endif
