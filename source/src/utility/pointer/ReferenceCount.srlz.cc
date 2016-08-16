// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/pointer/ReferenceCount.srlz.hh
/// @brief  Faux serialization routines for the ReferenceCount base class.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

// Unit headers
#include <utility/pointer/ReferenceCount.srlz.hh>

// Package headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>

namespace utility {
namespace pointer {

template < class Archive >
void save( Archive &, ReferenceCount const & ) {}

/// @brief External deserialization function for class ReferenceCount
template < class Archive >
void load( Archive &, ReferenceCount & ) {}

EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( ReferenceCount );

}
}

CEREAL_REGISTER_TYPE( utility::pointer::ReferenceCount )
CEREAL_REGISTER_DYNAMIC_INIT( utility_pointer_ReferenceCount )

#endif
