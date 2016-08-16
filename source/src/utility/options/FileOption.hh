// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/FileOption.hh
/// @brief  Program file option class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_FileOption_hh
#define INCLUDED_utility_options_FileOption_hh


// Unit headers
#include <utility/options/FileOption.fwd.hh>

// Package headers
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/keys/FileOptionKey.hh>

// Project headers
#include <utility/file/FileName.hh>

// ObjexxFCL headers
#include <ObjexxFCL/char.functions.hh>


namespace utility {
namespace options {


/// @brief Program file option class
class FileOption :
	public ScalarOption_T_< FileOptionKey, file::FileName >
{


private: // Types


	typedef  ScalarOption_T_< FileOptionKey, file::FileName >  Super;
	typedef  file::FileName  FileName;


public: // Types


	using Super::default_value;
	using Super::def;
	using Super::value;
	using Super::operator ();


public: // Creation


	/// @brief Default constructor
	inline
	FileOption()
	{}


	/// @brief Key + description constructor
	inline
	FileOption(
		FileOptionKey const & key_a,
		std::string const & description_a
	) :
		Super( key_a, description_a )
	{}


	/// @brief Clone this
	inline
	FileOption *
	clone() const
	{
		return new FileOption( *this );
	}


	/// @brief Destructor
	inline
	virtual
	~FileOption()
	{}


public: // Methods


	/// @brief Default value assignment
	inline
	FileOption &
	default_value( std::string const & value_a )
	{
		Super::default_value( FileName( value_a ) );
		return *this;
	}


	/// @brief Default value assignment
	inline
	FileOption &
	def( std::string const & value_a )
	{
		Super::default_value( FileName( value_a ) );
		return *this;
	}


	/// @brief Value assignment
	inline
	FileOption &
	value( std::string const & value_a )
	{
		Super::value( FileName( value_a ) );
		return *this;
	}


	/// @brief Value assignment
	inline
	FileOption &
	operator ()( std::string const & value_a )
	{
		Super::operator ()( FileName( value_a ) );
		return *this;
	}

	/// @brief function to allow atomatic conversion to the string.
	inline
	operator std::string() { return this->value(); }


	/// @brief function to allow atomatic conversion to the char *.
	typedef char const * CHAR_CONST_P;
	inline
	operator CHAR_CONST_P() { cached_value_=this->value(); return cached_value_.c_str(); }


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
		return "F";
	}


protected: // Methods


	/// @brief Value of a string
	inline
	Value
	value_of( std::string const & value_str ) const
	{
		return FileName( value_str );
	}

private: // data

	/// @brief temp storage for option value in string form. We need this to make type conversion to
	/// char * possible.
	std::string cached_value_;

}; // FileOption


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_FileOption_HH
