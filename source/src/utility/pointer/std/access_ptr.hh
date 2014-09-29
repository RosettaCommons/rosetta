// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/std/access_ptr.hh
/// @brief  Non-owning access smart pointer
/// @author Luki Goldschmidt <lugo@uw.edu>


#ifndef INCLUDED_utility_pointer_std_access_ptr_hh
#define INCLUDED_utility_pointer_std_access_ptr_hh

#ifdef PTR_STD

// Unit headers
#include <utility/pointer/owning_ptr.hh>

// Project headers
#include <utility/down_cast.hh>

// C++ headers
#include <cassert>
#include <iosfwd>
#include <memory>

namespace utility {
namespace pointer {

using std::weak_ptr;
using std::owner_less;

/// @brief Equality comparator
template< typename T, typename U >
inline
bool
equal( weak_ptr< T > const & a, weak_ptr< U > const & b )
{
	shared_ptr< T > as = a.lock();
	shared_ptr< U > bs = b.lock();
	return as && bs && as.get() == bs.get();
}

/// @brief Equality comparator
template< typename T, typename U >
inline
bool
equal( weak_ptr< T > & a, shared_ptr< U > const & bs )
{
	shared_ptr< T > as = a.lock();
	return as && as.get() == bs.get();
}

/// @brief Equality comparator
template< typename T, typename U >
inline
bool
equal( weak_ptr< T > & a,  U* const b )
{
	shared_ptr< T > as = a.lock();
	return as && as.get() == b;
}

} // namespace pointer
} // namespace utility

#endif // PTR_STD


#endif // INCLUDED_utility_pointer_std_access_ptr_HH
