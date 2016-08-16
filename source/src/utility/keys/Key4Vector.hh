// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/Key4Vector.hh
/// @brief  4-key meta-key
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note
///  @li Provides a meta-key from keys of the same type


#ifndef INCLUDED_utility_keys_Key4Vector_hh
#define INCLUDED_utility_keys_Key4Vector_hh


// Unit headers
#include <utility/keys/Key4Vector.fwd.hh>

// C++ headers
#include <utility/assert.hh>


namespace utility {
namespace keys {


/// @brief 4-key meta-key
template< typename K >
class Key4Vector
{


public: // Types


	typedef  K  Key;


public: // Creation


	/// @brief Default constructor
	/// @note  Only works if Key has default constructor
	inline
	Key4Vector()
	{}


	/// @brief Key constructor
	inline
	Key4Vector(
		Key const & key1_a,
		Key const & key2_a,
		Key const & key3_a,
		Key const & key4_a
	) :
		key1_( key1_a ),
		key2_( key2_a ),
		key3_( key3_a ),
		key4_( key4_a )
	{}


	/// @brief Destructor
	inline
	~Key4Vector()
	{}


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


	/// @brief Key 3
	inline
	Key const &
	key3() const
	{
		return key3_;
	}


	/// @brief Key 3
	inline
	Key &
	key3()
	{
		return key3_;
	}


	/// @brief Key 4
	inline
	Key const &
	key4() const
	{
		return key4_;
	}


	/// @brief Key 4
	inline
	Key &
	key4()
	{
		return key4_;
	}


public: // Indexers


	/// @brief Key4Vector[ i ] const: 0-based index
	inline
	Key const &
	operator []( int const i ) const
	{
	debug_assert( ( i >= 0 ) && ( i < 4 ) );
		return (
			( i == 0 ? key1_ :
			( i == 1 ? key2_ :
			( i == 2 ? key3_ : key4_ ) ) ) );
	}


	/// @brief Key4Vector[ i ]: 0-based index
	inline
	Key &
	operator []( int const i )
	{
	debug_assert( ( i >= 0 ) && ( i < 4 ) );
		return (
			( i == 0 ? key1_ :
			( i == 1 ? key2_ :
			( i == 2 ? key3_ : key4_ ) ) ) );
	}


	/// @brief Key4Vector( i ) const: 1-based index
	inline
	Key const &
	operator ()( int const i ) const
	{
	debug_assert( ( i > 0 ) && ( i <= 4 ) );
		return (
			( i == 1 ? key1_ :
			( i == 2 ? key2_ :
			( i == 3 ? key3_ : key4_ ) ) ) );
	}


	/// @brief Key4Vector( i ): 1-based index
	inline
	Key &
	operator ()( int const i )
	{
	debug_assert( ( i > 0 ) && ( i <= 4 ) );
		return (
			( i == 1 ? key1_ :
			( i == 2 ? key2_ :
			( i == 3 ? key3_ : key4_ ) ) ) );
	}


public: // Comparison


	/// @brief Key4Vector == Key4Vector
	friend
	inline
	bool
	operator ==( Key4Vector const & a, Key4Vector const & b )
	{
		return (
			( a.key1_ == b.key1_ ) &&
			( a.key2_ == b.key2_ ) &&
			( a.key3_ == b.key3_ ) &&
			( a.key4_ == b.key4_ ) );
	}


	/// @brief Key4Vector != Key4Vector
	friend
	inline
	bool
	operator !=( Key4Vector const & a, Key4Vector const & b )
	{
		return (
			( a.key1_ != b.key1_ ) ||
			( a.key2_ != b.key2_ ) ||
			( a.key3_ != b.key3_ ) ||
			( a.key4_ != b.key4_ ) );
	}


	/// @brief Key4Vector < Key4Vector
	/// @note  Lexicographic (full) ordering => Key4Vector is suitable for use as a map key or set element
	friend
	inline
	bool
	operator <( Key4Vector const & a, Key4Vector const & b )
	{
		return (
		 ( a.key1_ < b.key1_ ? true :
		 ( b.key1_ < a.key1_ ? false : // a.key1_ == b.key1_
		 ( a.key2_ < b.key2_ ? true :
		 ( b.key2_ < a.key2_ ? false : // a.key2_ == b.key2_
		 ( a.key3_ < b.key3_ ? true :
		 ( b.key3_ < a.key3_ ? false : // a.key3_ == b.key3_
		 ( a.key4_ < b.key4_ ) ) ) ) ) ) ) );
	}


private: // Fields


	/// @brief Keys
	Key key1_;
	Key key2_;
	Key key3_;
	Key key4_;


}; // Key4Vector


// Friend function namespace declarations


/// @brief Key4Vector == Key4Vector
template< typename K >
bool
operator ==( Key4Vector< K > const & a, Key4Vector< K > const & b );


/// @brief Key4Vector != Key4Vector
template< typename K >
bool
operator !=( Key4Vector< K > const & a, Key4Vector< K > const & b );


/// @brief Key4Vector < Key4Vector
template< typename K >
bool
operator <( Key4Vector< K > const & a, Key4Vector< K > const & b );


} // namespace keys
} // namespace utility


#endif // INCLUDED_utility_keys_Key4Vector_HH
