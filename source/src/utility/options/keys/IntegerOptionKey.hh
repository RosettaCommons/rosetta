// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/keys/IntegerOptionKey.hh
/// @brief  Automatic hidden index key for integer options
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_keys_IntegerOptionKey_hh
#define INCLUDED_utility_options_keys_IntegerOptionKey_hh


// Unit headers
#include <utility/options/keys/IntegerOptionKey.fwd.hh>

// Package headers
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/IntegerOption.fwd.hh>

// Project headers
#include <utility/keys/SmallKeyVector.fwd.hh>


namespace utility {
namespace options {


/// @brief Automatic hidden index key for integer options
class IntegerOptionKey :
	public ScalarOptionKey
{


private: // Types


	typedef  ScalarOptionKey  Super;


private: // Friends


	friend class utility::keys::SmallKeyVector< IntegerOptionKey, IntegerOption >;


public: // Creation


	/// @brief Default constructor
	inline
	IntegerOptionKey()
	{}


	/// @brief Copy + identifier constructor
	inline
	IntegerOptionKey(
		IntegerOptionKey const & key,
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
	IntegerOptionKey( Key const & key ) :
		Super( key )
	{}


	/// @brief Key + identifier constructor
	inline
	IntegerOptionKey(
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
	IntegerOptionKey(
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
	IntegerOptionKey *
	clone() const
	{
		return new IntegerOptionKey( *this );
	}


	/// @brief Destructor
	inline
	virtual
	~IntegerOptionKey() {}


public: // Assignment


	/// @brief Key assignment
	inline
	IntegerOptionKey &
	operator =( Key const & key )
	{
		assign_Key( key );
		return *this;
	}


}; // IntegerOptionKey


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_keys_IntegerOptionKey_HH
