// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/ReferenceCount.hh
/// @brief  Base class for reference-counted single inheritance polymorphic classes
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


// Unit headers
#include <utility/pointer/ReferenceCount.hh>


namespace utility {
namespace pointer {


/// @brief ReferenceCount static member definitions

#ifdef MULTI_THREADED
long const ReferenceCount::max_count_ = 123456789;
#else
ReferenceCount::Size const ReferenceCount::max_count_ = static_cast< ReferenceCount::Size >( -1 );
#endif

} // namespace pointer
} // namespace utility
