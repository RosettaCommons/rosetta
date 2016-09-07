// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/IntegerOption.hh
/// @brief  Program integer option class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_IntegerOption_hh
#define INCLUDED_utility_options_IntegerOption_hh


// Unit headers
#include <utility/options/IntegerOption.fwd.hh>

// Package headers
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/mpi_stderr.hh>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <cstdlib>
#include <iostream>


namespace utility {
namespace options {

/// @brief Program integer option class
class IntegerOption :
	public ScalarOption_T_< IntegerOptionKey, int >
{


private: // Types


	typedef  ScalarOption_T_< IntegerOptionKey, int >  Super;


public: // Creation


	/// @brief Default constructor
	inline
	IntegerOption()
	{}


	/// @brief Key + description constructor
	inline
	IntegerOption(
		IntegerOptionKey const & key_a,
		std::string const & description_a
	) :
		Super( key_a, description_a )
	{}


	/// @brief Clone this
	inline
	IntegerOption *
	clone() const override
	{
		return new IntegerOption( *this );
	}


	/// @brief Destructor
	inline
	
	~IntegerOption() override
	= default;


public: // Properties


	/// @brief Is a string readable as this option's value type?
	inline
	bool
	is_value( std::string const & value_str ) const override
	{
		return ObjexxFCL::is_int( value_str );
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
		return "I";
	}


protected: // Methods


	/// @brief Value of a string
	inline
	Value
	value_of( std::string const & value_str ) const override
	{
		if ( ( value_str.empty() ) || ( ! ObjexxFCL::is_int( value_str ) ) ) {
			mpi_safe_std_err("ERROR: Illegal value for integer option -" +id()+ " specified: " + value_str);
			std::exit( EXIT_FAILURE );
		}
		return ObjexxFCL::int_of( value_str );
	}


}; // IntegerOption


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_IntegerOption_HH
