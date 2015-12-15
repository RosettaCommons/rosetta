// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/ReferenceCount.srlz.hh
/// @brief  Faux serialization routines for the ReferenceCount base class.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

#ifndef INCLUDED_utility_pointer_ReferenceCount_SRLZ_HH
#define INCLUDED_utility_pointer_ReferenceCount_SRLZ_HH

#include <utility/pointer/ReferenceCount.hh>

// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>

namespace utility {
namespace pointer {

/// @brief External serialization function for class ReferenceCount to trick classes
/// that contain ReferenceCountOPs into serializing them.  Of course, no one is going
/// to instantiate a ReferenceCount and put it into a ReferenceCountOP, but without a
/// pure virtual function in ReferenceCount, there's no way to tell the Cereal library
/// not to worry about that contingency.
template < class Archive >
void save( Archive & arc, ReferenceCount const & rc );

/// @brief External deserialization function for class ReferenceCount
template < class Archive >
void load( Archive & arc, ReferenceCount & rc );

}
}

CEREAL_FORCE_DYNAMIC_INIT( utility_pointer_ReferenceCount )

#endif
#endif
