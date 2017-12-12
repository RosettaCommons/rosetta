// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/Key4Tuple.hh
/// @brief  4-tuple meta-key
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
///
/// @note
///  @li Provides a meta-key from keys of different types


#ifndef INCLUDED_utility_keys_Key4Tuple_hh
#define INCLUDED_utility_keys_Key4Tuple_hh


// Unit headers
#include <utility/keys/Key4Tuple.fwd.hh>


namespace utility {
namespace keys {


/// @brief 4-tuple meta-key
template<
	typename K1,
	typename K2,
	typename K3,
	typename K4
>
class Key4Tuple
{


public: // Types


	typedef  K1  Key1;
	typedef  K2  Key2;
	typedef  K3  Key3;
	typedef  K4  Key4;


public: // Creation


	/// @brief Default constructor
	/// @note  Only works if Keys have default constructors
	inline
	Key4Tuple()
	= default;


	/// @brief Key constructor
	inline
	Key4Tuple(
		Key1 const & key1_a,
		Key2 const & key2_a,
		Key3 const & key3_a,
		Key4 const & key4_a
	) :
		key1_( key1_a ),
		key2_( key2_a ),
		key3_( key3_a ),
		key4_( key4_a )
	{}


	/// @brief Destructor
	inline
	~Key4Tuple()
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


	/// @brief Key 3
	inline
	Key3 const &
	key3() const
	{
		return key3_;
	}


	/// @brief Key 3
	inline
	Key3 &
	key3()
	{
		return key3_;
	}


	/// @brief Key 4
	inline
	Key4 const &
	key4() const
	{
		return key4_;
	}


	/// @brief Key 4
	inline
	Key4 &
	key4()
	{
		return key4_;
	}


public: // Comparison


	/// @brief Key4Tuple == Key4Tuple
	friend
	inline
	bool
	operator ==( Key4Tuple const & a, Key4Tuple const & b )
	{
		return (
			( a.key1_ == b.key1_ ) &&
			( a.key2_ == b.key2_ ) &&
			( a.key3_ == b.key3_ ) &&
			( a.key4_ == b.key4_ ) );
	}


	/// @brief Key4Tuple != Key4Tuple
	friend
	inline
	bool
	operator !=( Key4Tuple const & a, Key4Tuple const & b )
	{
		return (
			( a.key1_ != b.key1_ ) ||
			( a.key2_ != b.key2_ ) ||
			( a.key3_ != b.key3_ ) ||
			( a.key4_ != b.key4_ ) );
	}


	/// @brief Key4Tuple < Key4Tuple
	/// @note  Lexicographic (full) ordering => Key4Tuple is suitable for use as a map key or set element
	friend
	inline
	bool
	operator <( Key4Tuple const & a, Key4Tuple const & b )
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
	Key1 key1_;
	Key2 key2_;
	Key3 key3_;
	Key4 key4_;


}; // Key4Tuple


// Friend function namespace declarations


/// @brief Key4Tuple == Key4Tuple
template<
	typename K1,
	typename K2,
	typename K3,
	typename K4
>
bool
operator ==( Key4Tuple< K1, K2, K3, K4 > const & a, Key4Tuple< K1, K2, K3, K4 > const & b );


/// @brief Key4Tuple != Key4Tuple
template<
	typename K1,
	typename K2,
	typename K3,
	typename K4
>
bool
operator !=( Key4Tuple< K1, K2, K3, K4 > const & a, Key4Tuple< K1, K2, K3, K4 > const & b );


/// @brief Key4Tuple < Key4Tuple
template<
	typename K1,
	typename K2,
	typename K3,
	typename K4
>
bool
operator <( Key4Tuple< K1, K2, K3, K4 > const & a, Key4Tuple< K1, K2, K3, K4 > const & b );


} // namespace keys
} // namespace utility


#endif // INCLUDED_utility_keys_Key4Tuple_HH
