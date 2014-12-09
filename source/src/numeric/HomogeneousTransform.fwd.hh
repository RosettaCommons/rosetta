// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/HomogeneousTransform.fwd.hh
/// @brief  Fast coordinate frame container
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
///
/// @remarks
///  @li Inline, loop-free functions for speed
///  @li Non-virtual destructor for speed: Not set up for use as a base class
///  @li Represents 4x4 homogenous matrix as a 4x3 table with the last row
///      implicitly represented as [ 0, 0, 0, 1 ]

#ifndef INCLUDED_numeric_HomogeneousTransform_fwd_hh
#define INCLUDED_numeric_HomogeneousTransform_fwd_hh

namespace numeric {

template  < class T >
class HomogeneousTransform;

}

#endif
