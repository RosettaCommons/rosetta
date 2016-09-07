// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/Key2Vector.hh
/// @brief  2-key meta-key
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note
///  @li Provides a meta-key from keys of the same type


#ifndef INCLUDED_utility_keys_Key2Vector_hh
#define INCLUDED_utility_keys_Key2Vector_hh


// Unit headers
#include <utility/keys/Key2Vector.fwd.hh>

// C++ headers
#include <utility/assert.hh>


namespace utility {
namespace keys {


/// @brief 2-key meta-key
template< typename K >
class Key2Vector
{


public: // Types


	typedef  K  Key;


public: // Creation


	/// @brief Default constructor
	/// @note  Only works if Key has default constructor
	inline
	Key2Vector()
	= default;


	/// @brief Key constructor
	inline
	Key2Vector(
		Key const & key1_a,
		Key const & key2_a
	) :
		key1_( key1_a ),
		key2_( key2_a )
	{}


	/// @brief Destructor
	inline
	~Key2Vector()
	= default;


public: // Properties


	/// @brief Key 1
	inline
	Key const &
	key1() const
	{
		return key1_;
	}


	/// @brief Key 1
	inline
	Key &
	key1()
	{
		return key1_;
	}


	/// @brief Key 2
	inline
	Key const &
	key2() const
	{
		return key2_;
	}


	/// @brief Key 2
	inline
	Key &
	key2()
	{
		return key2_;
	}


public: // Indexers


	/// @brief Key2Vector[ i ] const: 0-based index
	inline
	Key const &
	operator []( int const i ) const
	{
	debug_assert( ( i >= 0 ) && ( i < 2 ) );
		return ( i == 0 ? key1_ : key2_ );
	}


	/// @brief Key2Vector[ i ]: 0-based index
	inline
	Key &
	operator []( int const i )
	{
	debug_assert( ( i >= 0 ) && ( i < 2 ) );
		return ( i == 0 ? key1_ : key2_ );
	}


	/// @brief Key2Vector( i ) const: 1-based index
	inline
	Key const &
	operator ()( int const i ) const
	{
	debug_assert( ( i > 0 ) && ( i <= 2 ) );
		return ( i == 1 ? key1_ : key2_ );
	}


	/// @brief Key2Vector( i ): 1-based index
	inline
	Key &
	operator ()( int const i )
	{
	debug_assert( ( i > 0 ) && ( i <= 2 ) );
		return ( i == 1 ? key1_ : key2_ );
	}


public: // Comparison


	/// @brief Key2Vector == Key2Vector
	friend
	inline
	bool
	operator ==( Key2Vector const & a, Key2Vector const & b )
	{
		return (
			( a.key1_ == b.key1_ ) &&
			( a.key2_ == b.key2_ ) );
	}


	/// @brief Key2Vector != Key2Vector
	friend
	inline
	bool
	operator !=( Key2Vector const & a, Key2Vector const & b )
	{
		return (
			( a.key1_ != b.key1_ ) ||
			( a.key2_ != b.key2_ ) );
	}


	/// @brief Key2Vector < Key2Vector
	/// @note  Lexicographic (full) ordering => Key2Vector is suitable for use as a map key or set element
	friend
	inline
	bool
	operator <( Key2Vector const & a, Key2Vector const & b )
	{
		return (
			( a.key1_ < b.key1_ ? true :
			( b.key1_ < a.key1_ ? false : // a.key1_ == b.key1_
			( a.key2_ < b.key2_ ) ) ) );
	}


private: // Fields


	/// @brief Keys
	Key key1_;
	Key key2_;


}; // Key2Vector


// Friend function namespace declarations


/// @brief Key2Vector == Key2Vector
template< typename K >
bool
operator ==( Key2Vector< K > const & a, Key2Vector< K > const & b );


/// @brief Key2Vector != Key2Vector
template< typename K >
bool
operator !=( Key2Vector< K > const & a, Key2Vector< K > const & b );


/// @brief Key2Vector < Key2Vector
template< typename K >
bool
operator <( Key2Vector< K > const & a, Key2Vector< K > const & b );


} // namespace keys
} // namespace utility


#endif // INCLUDED_utility_keys_Key2Vector_HH
