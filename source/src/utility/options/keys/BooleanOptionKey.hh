// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/keys/BooleanOptionKey.hh
/// @brief  Automatic hidden index key for boolean options
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_keys_BooleanOptionKey_hh
#define INCLUDED_utility_options_keys_BooleanOptionKey_hh


// Unit headers
#include <utility/options/keys/BooleanOptionKey.fwd.hh>

// Package headers
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/BooleanOption.fwd.hh>

// Project headers
#include <utility/keys/SmallKeyVector.fwd.hh>


namespace utility {
namespace options {


/// @brief Automatic hidden index key for boolean options
class BooleanOptionKey :
	public ScalarOptionKey
{


private: // Types


	typedef  ScalarOptionKey  Super;


private: // Friends


	friend class utility::keys::SmallKeyVector< BooleanOptionKey, BooleanOption >;


public: // Creation


	/// @brief Default constructor
	inline
	BooleanOptionKey()
	{}


	/// @brief Copy + identifier constructor
	inline
	BooleanOptionKey(
		BooleanOptionKey const & key,
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
	BooleanOptionKey( Key const & key ) :
		Super( key )
	{}


	/// @brief Key + identifier constructor
	inline
	BooleanOptionKey(
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
	BooleanOptionKey(
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
	BooleanOptionKey *
	clone() const
	{
		return new BooleanOptionKey( *this );
	}


	/// @brief Destructor
	inline
	virtual
	~BooleanOptionKey() {}


public: // Assignment


	/// @brief Key assignment
	inline
	BooleanOptionKey &
	operator =( Key const & key )
	{
		assign_Key( key );
		return *this;
	}


}; // BooleanOptionKey


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_keys_BooleanOptionKey_HH
