// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/down_cast.hh
/// @brief  Fast polymorphic down-casting functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note A fast polymorphic down-cast when the cast is known to be valid
/// @note The cast validity is assert-checked in debug builds


#ifndef INCLUDED_utility_down_cast_hh
#define INCLUDED_utility_down_cast_hh


// C++ headers
#include <utility/assert.hh>


namespace utility {


/// @brief Meta-programming classes to provide the pointer type for down_cast
template< typename T > struct RawType        { typedef  T *  Pointer; };
template< typename T > struct RawType< T & > { typedef  T *  Pointer; };
template< typename T > struct RawType< T * > { typedef  T *  Pointer; };


/// @brief Fast assert-checked polymorphic down-cast: reference argument
///
/// @note Usage: down_cast< Type & > where Type can be const-qualified
/// @note For down-casting when you know the cast is valid
/// @note Can't use for hierarchies with virtual base classes
/// @note Assert intentionally won't compile if a virtual base class is present
template< class Target, class Source >
inline
Target
down_cast( Source & s )
{
	debug_assert( dynamic_cast< typename RawType< Target >::Pointer >( &s ) == &s );
	return static_cast< Target >( s );
}


/// @brief Fast assert-checked polymorphic down-cast: pointer argument
///
/// @note Usage: down_cast< Type * > where Type can be const-qualified
/// @note For down-casting when you know the cast is valid
/// @note Can't use for hierarchies with virtual base classes
/// @note Assert intentionally won't compile if a virtual base class is present
template< class Target, class Source >
inline
Target
down_cast( Source * p )
{
	debug_assert( dynamic_cast< typename RawType< Target >::Pointer >( p ) == p );
	return static_cast< Target >( p );
}


} // namespace utility


#endif // INCLUDED_utility_down_cast_HH
