// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/pointer/VirtualBase.srlz.hh
/// @brief  Faux serialization routines for the VirtualBase class.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifdef SERIALIZATION

#ifndef INCLUDED_utility_VirtualBase_SRLZ_HH
#define INCLUDED_utility_VirtualBase_SRLZ_HH

#include <utility/VirtualBase.hh>

// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>

namespace utility {

/// @brief External serialization function for class VirtualBase to trick classes
/// that contain VirtualBaseOPs into serializing them.  Of course, no one is going
/// to instantiate a VirtualBase and put it into a VirtualBaseOP, but without a
/// pure virtual function in VirtualBase, there's no way to tell the Cereal library
/// not to worry about that contingency.
template < class Archive >
void save( Archive & arc, VirtualBase const & rc );

/// @brief External deserialization function for class VirtualBase
template < class Archive >
void load( Archive & arc, VirtualBase & rc );

} // namespace utility

CEREAL_FORCE_DYNAMIC_INIT( utility_VirtualBase )

#endif // INCLUDED_utility_VirtualBase_SRLZ_HH
#endif // SERIALIZATION
