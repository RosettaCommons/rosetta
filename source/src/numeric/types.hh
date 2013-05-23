// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   numeric/types.hh
/// @brief  rosetta project type declarations. Should be kept updated with
///         core/types.hh. This exists because numeric can't depend on
///         files in core.
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James Thompson


#ifndef INCLUDED_numeric_types_hh
#define INCLUDED_numeric_types_hh


// Numeric headers
#include <numeric/xyzVector.fwd.hh>

// Platform headers
#include <platform/types.hh> // ssize_t

// C++ headers
#include <cstddef> // std::size_t


namespace numeric {


// Floating point precision control scalar
#ifdef ROSETTA_FLOAT // Real == float
typedef  float  Real;
#else // Real == double
typedef  double  Real;
#endif

typedef platform::Size Size;
typedef platform::SSize SSize;

} // namespace numeric


#endif // INCLUDED_numeric_types_HH
