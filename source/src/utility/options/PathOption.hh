// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/PathOption.hh
/// @brief  Program path option class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_PathOption_hh
#define INCLUDED_utility_options_PathOption_hh


// Unit headers
#include <utility/options/PathOption.fwd.hh>

// Package headers
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/keys/PathOptionKey.hh>

// Project headers
#include <utility/file/PathName.hh>

// ObjexxFCL headers
#include <ObjexxFCL/char.functions.hh>


namespace utility {
namespace options {


/// @brief Program path option class
class PathOption :
	public ScalarOption_T_< PathOptionKey, file::PathName >
{


private: // Types


	typedef  ScalarOption_T_< PathOptionKey, file::PathName >  Super;
	typedef  file::PathName  PathName;


public: // Types


	using Super::default_value;
	using Super::def;
	using Super::value;
	using Super::operator ();


public: // Creation


	/// @brief Default constructor
	inline
	PathOption()
	{}


	/// @brief Key + description constructor
	inline
	PathOption(
		PathOptionKey const & key_a,
		std::string const & description_a
	) :
		Super( key_a, description_a )
	{}


	/// @brief Clone this
	inline
	PathOption *
	clone() const
	{
		return new PathOption( *this );
	}


	/// @brief Destructor
	inline
	virtual
	~PathOption()
	{}


public: // Methods


	/// @brief Default value assignment
	inline
	PathOption &
	default_value( std::string const & value_a )
	{
		Super::default_value( PathName( value_a ) );
		return *this;
	}


	/// @brief Default value assignment
	inline
	PathOption &
	def( std::string const & value_a )
	{
		Super::default_value( PathName( value_a ) );
		return *this;
	}


	/// @brief Value assignment
	inline
	PathOption &
	value( std::string const & value_a )
	{
		Super::value( PathName( value_a ) );
		return *this;
	}


	/// @brief Value assignment
	inline
	PathOption &
	operator ()( std::string const & value_a )
	{
		Super::operator ()( PathName( value_a ) );
		return *this;
	}


public: // Properties


	/// @brief Is a string readable as this option's value type?
	inline
	bool
	is_value( std::string const & ) const
	{
		return true;
	}


	/// @brief Is a string readable as this option's value type and a legal command line value?
	inline
	bool
	is_cl_value( std::string const & value_str ) const
	{
		return ( ( value_str.empty() ) || ( ! ObjexxFCL::is_any_of( value_str[ 0 ], "-@" ) ) );
	}


	/// @brief Option type code string representation
	inline
	std::string
	type_string() const
	{
		return "P";
	}


protected: // Methods


	/// @brief Value of a string
	inline
	Value
	value_of( std::string const & value_str ) const
	{
		return PathName( value_str );
	}


}; // PathOption


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_PathOption_HH
