// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/keys/ScalarOptionKey.hh
/// @brief  Abstract automatic hidden index key for scalar-valued options
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_keys_ScalarOptionKey_hh
#define INCLUDED_utility_options_keys_ScalarOptionKey_hh


// Unit headers
#include <utility/options/keys/ScalarOptionKey.fwd.hh>

// Package headers
#include <utility/options/keys/OptionKey.hh>


namespace utility {
namespace options {


/// @brief Abstract automatic hidden index key for scalar-valued options
class ScalarOptionKey :
	public OptionKey
{


private: // Types


	typedef  OptionKey  Super;


protected: // Creation


	/// @brief Default constructor
	inline
	ScalarOptionKey()
	{}


	/// @brief Copy constructor
	inline
	ScalarOptionKey( ScalarOptionKey const & key ) :
		Super( key )
	{}


	/// @brief Copy + identifier constructor
	inline
	ScalarOptionKey(
		ScalarOptionKey const & key,
		std::string const & id_a,
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( key, id_a, identifier_a, code_a )
	{}


	/// @brief Key constructor
	inline
	explicit
	ScalarOptionKey( Key const & key ) :
		Super( key )
	{}


	/// @brief Key + identifier constructor
	inline
	ScalarOptionKey(
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
	ScalarOptionKey(
		std::string const & id_a,
		std::string const & identifier_a = std::string(),
		std::string const & code_a = std::string()
	) :
		Super( id_a, identifier_a, code_a )
	{}


public: // Creation


	/// @brief Clone this
	virtual
	ScalarOptionKey *
	clone() const = 0;


	/// @brief Destructor
	inline
	virtual
	~ScalarOptionKey()
	{}


public: // Assignment


	/// @brief Key assignment
	inline
	ScalarOptionKey &
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
		return true;
	}


	/// @brief Vector option key?
	inline
	bool
	vector() const
	{
		return false;
	}


}; // ScalarOptionKey


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_keys_ScalarOptionKey_HH
