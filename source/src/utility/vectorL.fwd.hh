// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/vectorL.fwd.hh
/// @brief  utility::vectorL forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_vectorL_fwd_hh
#define INCLUDED_utility_vectorL_fwd_hh


// Platform headers
#include <platform/types.hh>


// std::allocator Declaration
#ifdef UNUSUAL_ALLOCATOR_DECLARATION

// C++ headers
#include <vector>
#include <memory>

#else

// Faster but not 100% portable
namespace std { template< typename > class allocator; }

#endif // UNUSUAL_ALLOCATOR_DECLARATION


namespace utility {


// Forward
template< ssize_t, typename T, typename A = std::allocator< T > > class vectorL;


} // namespace utility


#endif // INCLUDED_utility_vectorL_FWD_HH
