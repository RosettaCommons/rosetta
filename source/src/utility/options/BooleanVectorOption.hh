// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/BooleanVectorOption.hh
/// @brief  Program boolean vector option class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_BooleanVectorOption_hh
#define INCLUDED_utility_options_BooleanVectorOption_hh


// Unit headers
#include <utility/options/BooleanVectorOption.fwd.hh>

// Package headers
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>

// C++ headers
#include <cstdlib>
#include <iostream>


namespace utility {
namespace options {


/// @brief Program boolean vector option class
class BooleanVectorOption :
	public VectorOption_T_< BooleanVectorOptionKey, bool >
{


private: // Types


	typedef  VectorOption_T_< BooleanVectorOptionKey, bool >  Super;


public: // Creation


	/// @brief Default constructor
	inline
	BooleanVectorOption()
	{}


	/// @brief Key + description constructor
	inline
	BooleanVectorOption(
		BooleanVectorOptionKey const & key_a,
		std::string const & description_a
	) :
		Super( key_a, description_a )
	{}


	/// @brief Clone this
	inline
	BooleanVectorOption *
	clone() const override
	{
		return new BooleanVectorOption( *this );
	}


	/// @brief Destructor
	inline

	~BooleanVectorOption() override
	= default;


public: // Properties


	/// @brief Is a string readable as this option's value type?
	inline
	bool
	is_value( std::string const & value_str ) const override
	{
		return ( ( is_true_value( value_str ) ) || ( is_false_value( value_str ) ) );
	}


	/// @brief Is a string readable as this option's value type and a legal command line value?
	inline
	bool
	is_cl_value( std::string const & value_str ) const override
	{
		return is_value( value_str );
	}


	/// @brief Option type code string representation
	inline
	std::string
	type_string() const override
	{
		return "(B" + size_constraint_string() + ')';
	}


	/// @brief Legal value string representation
	inline
	std::string
	legal_string() const override
	{
		return std::string();
	}


protected: // Methods


	/// @brief Value of a string
	inline
	Value
	value_of( std::string const & value_str ) const override
	{
		if ( is_true_value( value_str ) ) {
			return true;
		} else if ( is_false_value( value_str ) ) {
			return false;
		} else { // Illegal
			std::cerr << "ERROR: Illegal value for boolean option -" << id()
				<< " specified: " << value_str << std::endl;
			std::exit( EXIT_FAILURE );
			return false; // Keep compiler happy
		}
	}


	/// @brief String accepted as a true value?
	inline
	bool
	is_true_value( std::string const & value_str ) const
	{
		return (
			( value_str.empty() ) ||
			( value_str == "true" ) ||
			( value_str == "True" ) ||
			( value_str == "TRUE" ) ||
			( value_str == "t" ) ||
			( value_str == "T" ) ||
			( value_str == "1" ) ||
			( value_str == "on" ) ||
			( value_str == "On" ) ||
			( value_str == "ON" ) ||
			( value_str == "yes" ) ||
			( value_str == "Yes" ) ||
			( value_str == "YES" ) );
	}


	/// @brief String accepted as a true value?
	inline
	bool
	is_false_value( std::string const & value_str ) const
	{
		return (
			( value_str == "false" ) ||
			( value_str == "False" ) ||
			( value_str == "FALSE" ) ||
			( value_str == "f" ) ||
			( value_str == "F" ) ||
			( value_str == "0" ) ||
			( value_str == "off" ) ||
			( value_str == "Off" ) ||
			( value_str == "OFF" ) ||
			( value_str == "no" ) ||
			( value_str == "No" ) ||
			( value_str == "NO" ) );
	}


}; // BooleanVectorOption


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_BooleanVectorOption_HH
