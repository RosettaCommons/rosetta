// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/keys/AnyVectorOptionKey.hh
/// @brief  Automatic hidden index key for any vector-valued options
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_keys_AnyVectorOptionKey_hh
#define INCLUDED_utility_options_keys_AnyVectorOptionKey_hh


// Unit headers
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>

// Package headers
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/AnyOption.fwd.hh>

// Project headers
#include <utility/keys/SmallKeyVector.fwd.hh>


namespace utility {
namespace options {


/// @brief Automatic hidden index key for any vector-valued options
class AnyVectorOptionKey :
	public VectorOptionKey
{


private: // Types


	typedef  VectorOptionKey  Super;


private: // Friends


#if !(defined _MSC_VER) || (defined __INTEL_COMPILER) // Visual C++ 2005 bug work-around
	template< typename K, typename T > friend class utility::keys::SmallKeyVector;
#endif


public: // Creation


	/// @brief Default constructor
	inline
	AnyVectorOptionKey()
	{}


	/// @brief Copy + identifier constructor
	inline
	AnyVectorOptionKey(
		AnyVectorOptionKey const & key,
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
	AnyVectorOptionKey( Key const & key ) :
		Super( key )
	{}


	/// @brief Key + identifier constructor
	inline
	AnyVectorOptionKey(
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
	AnyVectorOptionKey(
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
	AnyVectorOptionKey *
	clone() const
	{
		return new AnyVectorOptionKey( *this );
	}


	/// @brief Destructor
	inline
	virtual
	~AnyVectorOptionKey()
	{}


public: // Assignment


	/// @brief Key assignment
	inline
	AnyVectorOptionKey &
	operator =( Key const & key )
	{
		assign_Key( key );
		return *this;
	}


}; // AnyVectorOptionKey


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_keys_AnyVectorOptionKey_HH
