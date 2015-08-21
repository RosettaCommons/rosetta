// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/options/BooleanOption.hh
/// @brief  Program boolean option class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_BooleanOption_hh
#define INCLUDED_utility_options_BooleanOption_hh


// Unit headers
#include <utility/options/BooleanOption.fwd.hh>

// Package headers
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/keys/BooleanOptionKey.hh>

#include <utility/string_util.hh>

// C++ headers
#include <cstdlib>
#include <iostream>


namespace utility {
namespace options {


/// @brief Program boolean option class
class BooleanOption :
	public ScalarOption_T_< BooleanOptionKey, bool >
{


private: // Types


	typedef  ScalarOption_T_< BooleanOptionKey, bool >  Super;


public: // Creation


	/// @brief Default constructor
	inline
	BooleanOption()
	{}


	/// @brief Key + description constructor
	inline
	BooleanOption(
		BooleanOptionKey const & key_a,
		std::string const & description_a
	) :
		Super( key_a, description_a )
	{
		default_value( false ); // This is usually what you want
	}


	/// @brief Key + description no-default named constructor
	inline
	static
	BooleanOption
	NoDefault(
		BooleanOptionKey const & key_a,
		std::string const & description_a
	)
	{
		return BooleanOption( key_a, description_a, false );
	}


	/// @brief Clone this
	inline
	BooleanOption *
	clone() const
	{
		return new BooleanOption( *this );
	}


	/// @brief Destructor
	inline
	virtual
	~BooleanOption()
	{}


public: // Properties


	/// @brief Is a string readable as this option's value type?
	inline
	bool
	is_value( std::string const & value_str ) const
	{
		return ( ( is_true_value( value_str ) ) || ( is_false_value( value_str ) ) );
	}


	/// @brief Is a string readable as this option's value type and a legal command line value?
	inline
	bool
	is_cl_value( std::string const & value_str ) const
	{
		return is_value( value_str );
	}


	/// @brief Option type code string representation
	inline
	std::string
	type_string() const
	{
		return "B";
	}


	/// @brief Legal value string representation
	inline
	std::string
	legal_string() const
	{
		return std::string();
	}


	/// @brief Default value string representation
	inline
	std::string
	default_string() const
	{
		return std::string(); // Don't show boolean defaults
	}


	/// @brief Value string representation
	inline
	std::string
	value_string() const
	{
		return ( ( active() ) && ( value() == false ) ? "false" : "" );
	}


	/// @brief =Value string representation
	inline
	std::string
	equals_string() const
	{
		return ( active() ? ( value() == false ? "=false" : "" ) : "=" );
	}


protected: // Methods


	/// @brief Value of a string
	inline
	Value
	value_of( std::string const & value_str ) const
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
		return ( value_str.empty() || utility::is_true_string( value_str ) );
	}


	/// @brief String accepted as a false value?
	inline
	bool
	is_false_value( std::string const & value_str ) const
	{
		return utility::is_false_string( value_str );
	}

private: // Creation


	/// @brief Key + description default control constructor
	inline
	BooleanOption(
		BooleanOptionKey const & key_a,
		std::string const & description_a,
		bool const set_default
	) :
		Super( key_a, description_a )
	{
		if ( set_default ) default_value( false );
	}


}; // BooleanOption


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_BooleanOption_HH
