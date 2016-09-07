// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/StringVectorOption.hh
/// @brief  Program string vector option class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_StringVectorOption_hh
#define INCLUDED_utility_options_StringVectorOption_hh


// Unit headers
#include <utility/options/StringVectorOption.fwd.hh>

// Package headers
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>

// ObjexxFCL headers
#include <ObjexxFCL/char.functions.hh>


namespace utility {
namespace options {


/// @brief Program string option class
class StringVectorOption :
	public VectorOption_T_< StringVectorOptionKey, std::string >
{


private: // Types


	typedef  VectorOption_T_< StringVectorOptionKey, std::string >  Super;


public: // Creation


	/// @brief Default constructor
	inline
	StringVectorOption()
	{}


	/// @brief Key + description constructor
	inline
	StringVectorOption(
		StringVectorOptionKey const & key_a,
		std::string const & description_a
	) :
		Super( key_a, description_a )
	{}


	/// @brief Clone this
	inline
	StringVectorOption *
	clone() const override
	{
		return new StringVectorOption( *this );
	}


	/// @brief Destructor
	inline

	~StringVectorOption() override
	= default;


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
		return "(S" + size_constraint_string() + ')';
	}


protected: // Methods


	/// @brief Value of a string
	inline
	Value
	value_of( std::string const & value_str ) const override
	{
		return value_str;
	}


}; // StringVectorOption


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_StringVectorOption_HH
