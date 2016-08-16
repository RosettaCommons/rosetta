// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/KeyCount.hh
/// @brief  Key counter functor
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note   Useful for holding the number of keys in a namespace-scoped UserKey collection


#ifndef INCLUDED_utility_keys_KeyCount_hh
#define INCLUDED_utility_keys_KeyCount_hh


// Unit headers
#include <utility/keys/KeyCount.fwd.hh>

// C++ headers
#include <utility/assert.hh>
#include <cstddef>


namespace utility {
namespace keys {


/// @brief Key counter functor
class KeyCount
{


public: // Types


	// STL/boost style
	typedef  std::size_type  size_type;

	// Project style
	typedef  std::size_type  Size;


public: // Creation


	/// @brief Count constructor
	inline
	explicit
	KeyCount( Size const count_a ) :
		count_( count_a )
	{}


	/// @brief Count + expected count constructor
	/// @note  Useful if a namespace constant is stored with the number of keys
	///        so that other namespace-scoped UserKeys can set contiguous indexes
	///        without global initialization order issues
	inline
	explicit
	KeyCount( Size const count_a, Size const expected_count ) :
		count_( count_a )
	{
	debug_assert( count_ == expected_count );
		if ( this ); // Silly if to suppress unused variable warnings
	}


public: // Properties


	/// @brief Count
	inline
	Size
	operator ()() const
	{
		return count_;
	}


private: // Fields


	/// @brief Count of keys
	Size count_;


}; // KeyCount


} // namespace keys
} // namespace utility


#endif // INCLUDED_utility_keys_KeyCount_HH
