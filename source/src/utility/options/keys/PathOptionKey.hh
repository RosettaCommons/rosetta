// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/keys/PathOptionKey.hh
/// @brief  Automatic hidden index key for path options
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_keys_PathOptionKey_hh
#define INCLUDED_utility_options_keys_PathOptionKey_hh


// Unit headers
#include <utility/options/keys/PathOptionKey.fwd.hh>

// Package headers
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/PathOption.fwd.hh>

// Project headers
#include <utility/keys/SmallKeyVector.fwd.hh>


namespace utility {
namespace options {


/// @brief Automatic hidden index key for path options
class PathOptionKey :
	public ScalarOptionKey
{


private: // Types


	typedef  ScalarOptionKey  Super;


private: // Friends


	friend class utility::keys::SmallKeyVector< PathOptionKey, PathOption >;


public: // Creation


	/// @brief Default constructor
	inline
	PathOptionKey()
	{}


	/// @brief Copy + identifier constructor
	inline
	PathOptionKey(
		PathOptionKey const & key,
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
	PathOptionKey( Key const & key ) :
		Super( key )
	{}


	/// @brief Key + identifier constructor
	inline
	PathOptionKey(
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
	PathOptionKey(
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
	PathOptionKey *
	clone() const
	{
		return new PathOptionKey( *this );
	}


	/// @brief Destructor
	inline
	virtual
	~PathOptionKey()
	{}


public: // Assignment


	/// @brief Key assignment
	inline
	PathOptionKey &
	operator =( Key const & key )
	{
		assign_Key( key );
		return *this;
	}


}; // PathOptionKey


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_keys_PathOptionKey_HH
