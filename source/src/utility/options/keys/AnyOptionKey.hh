// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/keys/AnyOptionKey.hh
/// @brief  Automatic hidden index key for any scalar-valued options
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_keys_AnyOptionKey_hh
#define INCLUDED_utility_options_keys_AnyOptionKey_hh


// Unit headers
#include <utility/options/keys/AnyOptionKey.fwd.hh>

// Package headers
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/AnyOption.fwd.hh>

// Project headers
#include <utility/keys/SmallKeyVector.fwd.hh>


namespace utility {
namespace options {


/// @brief Automatic hidden index key for any scalar-valued options
class AnyOptionKey :
	public ScalarOptionKey
{


private: // Types


	typedef  ScalarOptionKey  Super;


private: // Friends


#if !(defined _MSC_VER) || (defined __INTEL_COMPILER) // Visual C++ 2005 bug work-around
	template< typename K, typename T > friend class utility::keys::SmallKeyVector;
#endif


public: // Creation


	/// @brief Default constructor
	inline
	AnyOptionKey()
	{}


	/// @brief Copy + identifier constructor
	inline
	AnyOptionKey(
		AnyOptionKey const & key,
		std::string const & id_a,
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( key, id_a, identifier_a, code_a )
	{
		Lookup::add( *this ); // Add key to lookup map
	}


	/// @brief Key constructor
	inline
	explicit
	AnyOptionKey( Key const & key ) :
		Super( key )
	{}


	/// @brief Key + identifier constructor
	inline
	AnyOptionKey(
		Key const & key,
		std::string const & id_a,
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( key, id_a, identifier_a, code_a )
	{
		Lookup::add( *this ); // Add key to lookup map
	}


	/// @brief Identifier constructor
	inline
	explicit
	AnyOptionKey(
		std::string const & id_a,
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( id_a, identifier_a, code_a )
	{
		Lookup::add( *this ); // Add key to lookup map
	}


	/// @brief Clone this
	inline
	AnyOptionKey *
	clone() const
	{
		return new AnyOptionKey( *this );
	}


	/// @brief Destructor
	inline
	virtual
	~AnyOptionKey()
	{}


public: // Assignment


	/// @brief Key assignment
	inline
	AnyOptionKey &
	operator =( Key const & key )
	{
		assign_Key( key );
		return *this;
	}


}; // AnyOptionKey


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_keys_AnyOptionKey_HH
