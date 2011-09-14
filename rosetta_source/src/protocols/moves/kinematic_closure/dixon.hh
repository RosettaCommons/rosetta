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

#ifndef INCLUDE_src_protocols_moves_kinematic_closure_dixon_hh
#define INCLUDE_src_protocols_moves_kinematic_closure_dixon_hh

// Rosetta Headers
#include <core/types.hh>

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>

//Auto Headers
#include <utility/vector1.fwd.hh>


namespace protocols {
namespace moves {
namespace kinematic_closure {

void dixon(const utility::vector1<utility::vector1<core::Real> >& A,
		   const utility::vector1<utility::vector1<core::Real> >& B,
		   const utility::vector1<utility::vector1<core::Real> >& C,
		   const utility::vector1<utility::vector1<core::Real> >& D,
		   const utility::vector1<int>& order, utility::vector1<utility::vector1<core::Real> >& cos,
		   utility::vector1<utility::vector1<core::Real> >& sin, utility::vector1<utility::vector1<core::Real> >& tau,
		   int& nsol);

} // end namespace kinematic_closure
} // end namespace moves
} // end namespace protocols

#endif
