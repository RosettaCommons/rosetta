// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/Key2Tuple.hh
/// @brief  2-tuple meta-key
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note
///  @li Provides a meta-key from keys of different types


#ifndef INCLUDED_utility_keys_Key2Tuple_hh
#define INCLUDED_utility_keys_Key2Tuple_hh


// Unit headers
#include <utility/keys/Key2Tuple.fwd.hh>


namespace utility {
namespace keys {


/// @brief 2-tuple meta-key
template<
	typename K1,
	typename K2
>
class Key2Tuple
{


public: // Types


	typedef  K1  Key1;
	typedef  K2  Key2;


public: // Creation


	/// @brief Default constructor
	/// @note  Only works if Keys have default constructors
	inline
	Key2Tuple()
	= default;


	/// @brief Key constructor
	inline
	Key2Tuple(
		Key1 const & key1_a,
		Key2 const & key2_a
	) :
		key1_( key1_a ),
		key2_( key2_a )
	{}


	/// @brief Destructor
	inline
	~Key2Tuple()
	= default;


public: // Properties


	/// @brief Key 1
	inline
	Key1 const &
	key1() const
	{
		return key1_;
	}


	/// @brief Key 1
	inline
	Key1 &
	key1()
	{
		return key1_;
	}


	/// @brief Key 2
	inline
	Key2 const &
	key2() const
	{
		return key2_;
	}


	/// @brief Key 2
	inline
	Key2 &
	key2()
	{
		return key2_;
	}


public: // Comparison


	/// @brief Key2Tuple == Key2Tuple
	friend
	inline
	bool
	operator ==( Key2Tuple const & a, Key2Tuple const & b )
	{
		return (
			( a.key1_ == b.key1_ ) &&
			( a.key2_ == b.key2_ ) );
	}


	/// @brief Key2Tuple != Key2Tuple
	friend
	inline
	bool
	operator !=( Key2Tuple const & a, Key2Tuple const & b )
	{
		return (
			( a.key1_ != b.key1_ ) ||
			( a.key2_ != b.key2_ ) );
	}


	/// @brief Key2Tuple < Key2Tuple
	/// @note  Lexicographic (full) ordering => Key2Tuple is suitable for use as a map key or set element
	friend
	inline
	bool
	operator <( Key2Tuple const & a, Key2Tuple const & b )
	{
		return (
			( a.key1_ < b.key1_ ? true :
			( b.key1_ < a.key1_ ? false : // a.key1_ == b.key1_
			( a.key2_ < b.key2_ ) ) ) );
	}


private: // Fields


	/// @brief Keys
	Key1 key1_;
	Key2 key2_;


}; // Key2Tuple


// Friend function namespace declarations


/// @brief Key2Tuple == Key2Tuple
template<
	typename K1,
	typename K2
>
bool
operator ==( Key2Tuple< K1, K2 > const & a, Key2Tuple< K1, K2 > const & b );


/// @brief Key2Tuple != Key2Tuple
template<
	typename K1,
	typename K2
>
bool
operator !=( Key2Tuple< K1, K2 > const & a, Key2Tuple< K1, K2 > const & b );


/// @brief Key2Tuple < Key2Tuple
/// @note  Lexicographic (full) ordering => Key2Tuple is suitable for use as a map key or set element
template<
	typename K1,
	typename K2
>
bool
operator <( Key2Tuple< K1, K2 > const & a, Key2Tuple< K1, K2 > const & b );


} // namespace keys
} // namespace utility


#endif // INCLUDED_utility_keys_Key2Tuple_HH
