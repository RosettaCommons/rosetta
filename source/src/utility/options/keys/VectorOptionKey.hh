// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/keys/VectorOptionKey.hh
/// @brief  Abstract automatic hidden index key for vector-valued options
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_keys_VectorOptionKey_hh
#define INCLUDED_utility_options_keys_VectorOptionKey_hh


// Unit headers
#include <utility/options/keys/VectorOptionKey.fwd.hh>

// Package headers
#include <utility/options/keys/OptionKey.hh>


namespace utility {
namespace options {


/// @brief Abstract automatic hidden index key for vector-valued options
class VectorOptionKey :
	public OptionKey
{


private: // Types


	typedef  OptionKey  Super;


protected: // Creation


	/// @brief Default constructor
	inline
	VectorOptionKey()
	{}


	/// @brief Copy constructor
	inline
	VectorOptionKey( VectorOptionKey const & key ) :
		Super( key )
	{}


	/// @brief Copy + identifier constructor
	inline
	VectorOptionKey(
		VectorOptionKey const & key,
		std::string const & id_a,
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( key, id_a, identifier_a, code_a )
	{}


	/// @brief Key constructor
	inline
	explicit
	VectorOptionKey( Key const & key ) :
		Super( key )
	{}


	/// @brief Key + identifier constructor
	inline
	VectorOptionKey(
		Key const & key,
		std::string const & id_a,
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( key, id_a, identifier_a, code_a )
	{}


	/// @brief Identifier constructor
	inline
	explicit
	VectorOptionKey(
		std::string const & id_a,
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( id_a, identifier_a, code_a )
	{}


public: // Creation


	/// @brief Clone this
	virtual
	VectorOptionKey *
	clone() const = 0;


	/// @brief Destructor
	inline
	virtual
	~VectorOptionKey()
	{}


public: // Assignment


	/// @brief Key assignment
	inline
	VectorOptionKey &
	operator =( Key const & key )
	{
		assign_Key( key );
		return *this;
	}


public: // Properties


	/// @brief Scalar option key?
	inline
	bool
	scalar() const
	{
		return false;
	}


	/// @brief Vector option key?
	inline
	bool
	vector() const
	{
		return true;
	}


}; // VectorOptionKey


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_keys_VectorOptionKey_HH
