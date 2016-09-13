// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/StringOption.hh
/// @brief  Program string option class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_StringOption_hh
#define INCLUDED_utility_options_StringOption_hh


// Unit headers
#include <utility/options/StringOption.fwd.hh>

// Package headers
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/keys/StringOptionKey.hh>

// ObjexxFCL headers
#include <ObjexxFCL/char.functions.hh>


namespace utility {
namespace options {


/// @brief Program string option class
class StringOption :
	public ScalarOption_T_< StringOptionKey, std::string >
{


private: // Types


	typedef  ScalarOption_T_< StringOptionKey, std::string >  Super;


public: // Creation


	/// @brief Default constructor
	inline
	StringOption()
	{}


	/// @brief Key + description constructor
	inline
	StringOption(
		StringOptionKey const & key_a,
		std::string const & description_a
	) :
		Super( key_a, description_a )
	{}


	/// @brief Clone this
	inline
	StringOption *
	clone() const override
	{
		return new StringOption( *this );
	}


	/// @brief Destructor
	inline

	~StringOption() {}


public: // Properties


	/// @brief Is a string readable as this option's value type?
	inline
	bool
	is_value( std::string const & ) const override
	{
		return true;
	}


	/// @brief Is a string readable as this option's value type and a legal command line value?
	inline
	bool
	is_cl_value( std::string const & value_str ) const override
	{
		return ( ( value_str.empty() ) || ( ! ObjexxFCL::is_any_of( value_str[ 0 ], "-@" ) ) );
	}


	/// @brief Option type code string representation
	inline
	std::string
	type_string() const override
	{
		return "S";
	}


public: // Concatenation


	/// @brief StringOption + StringOption
	friend
	inline
	std::string
	operator +( StringOption const & option1, StringOption const & option2 )
	{
		return option1() + option2();
	}


	/// @brief StringOption + std::string
	friend
	inline
	std::string
	operator +( StringOption const & option, std::string const & s )
	{
		return option() + s;
	}


	/// @brief std::string + StringOption
	friend
	inline
	std::string
	operator +( std::string const & s, StringOption const & option )
	{
		return s + option();
	}


protected: // Methods


	/// @brief Value of a string
	inline
	Value
	value_of( std::string const & value_str ) const override
	{
		return value_str;
	}


}; // StringOption


// Friend function namespace declarations


/// @brief StringOption + StringOption
#ifndef __clang__
std::string
operator +( StringOption const & option1, StringOption const & option2 );
#endif

/// @brief StringOption + std::string
#ifndef __clang__
std::string
operator +( StringOption const & option, std::string const & s );
#endif

/// @brief std::string + StringOption
#ifndef __clang__
std::string
operator +( std::string const & s, StringOption const & option );
#endif

} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_StringOption_HH
