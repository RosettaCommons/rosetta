// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/keys/RealVectorOptionKey.hh
/// @brief  Automatic hidden index key for real options
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_keys_RealVectorOptionKey_hh
#define INCLUDED_utility_options_keys_RealVectorOptionKey_hh


// Unit headers
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>

// Package headers
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/RealVectorOption.fwd.hh>

// Project headers
#include <utility/keys/SmallKeyVector.fwd.hh>


namespace utility {
namespace options {


/// @brief Automatic hidden index key for real options
class RealVectorOptionKey :
	public VectorOptionKey
{


private: // Types


	typedef  VectorOptionKey  Super;


private: // Friends


	friend class utility::keys::SmallKeyVector< RealVectorOptionKey, RealVectorOption >;


public: // Creation


	/// @brief Default constructor
	inline
	RealVectorOptionKey()
	{}


	/// @brief Copy + identifier constructor
	inline
	RealVectorOptionKey(
		RealVectorOptionKey const & key,
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
	RealVectorOptionKey( Key const & key ) :
		Super( key )
	{}


	/// @brief Key + identifier constructor
	inline
	RealVectorOptionKey(
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
	RealVectorOptionKey(
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
	RealVectorOptionKey *
	clone() const
	{
		return new RealVectorOptionKey( *this );
	}


	/// @brief Destructor
	inline
	virtual
	~RealVectorOptionKey()
	{}


public: // Assignment


	/// @brief Key assignment
	inline
	RealVectorOptionKey &
	operator =( Key const & key )
	{
		assign_Key( key );
		return *this;
	}


}; // RealVectorOptionKey


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_keys_RealVectorOptionKey_HH
