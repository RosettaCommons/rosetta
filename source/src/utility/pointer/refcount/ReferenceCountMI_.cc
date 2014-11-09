// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/ReferenceCountMI_.cc
/// @brief  Base class for reference-counted multiple inheritance polymorphic classes
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)

#ifdef PTR_REFCOUNT

// Unit headers
#include <utility/pointer/refcount/ReferenceCountMI_.hh>


namespace utility {
namespace pointer {


/// @brief ReferenceCountMI_ static member definitions
ReferenceCountMI_::Size const ReferenceCountMI_::max_count_ = static_cast< ReferenceCountMI_::Size >( -1 );


} // namespace pointer
} // namespace utility

#endif // PTR_REFCOUNT
