// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/serialization/ObjexxFCL/DynamicIndexRange.srlz.hh
/// @brief  Serlialization routines for ObjexxFCL::DynamicIndexRanges; these functions have to live in namespace
///         ObjexxFCL, even though they are being written and compiled in the utility library.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_serialization_ObjexxFCL_DynamicIndexRange_srlz_HH
#define INCLUDED_utility_serialization_ObjexxFCL_DynamicIndexRange_srlz_HH

#ifdef SERIALIZATION

// Unit headers
#include <ObjexxFCL/DynamicIndexRange.hh>

namespace ObjexxFCL {

template < class Archive >
void save( Archive & archive, ObjexxFCL::DynamicIndexRange const & ir );

template < class Archive >
void load( Archive & archive, ObjexxFCL::DynamicIndexRange & ir );

}

#endif // SERIALIZATION

#endif // INCLUDED_utility_serialization_ObjexxFCL_DynamicIndexRange_srlz_HH
