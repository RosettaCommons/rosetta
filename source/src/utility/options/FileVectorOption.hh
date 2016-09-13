// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/FileVectorOption.hh
/// @brief  Program file vector option class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_FileVectorOption_hh
#define INCLUDED_utility_options_FileVectorOption_hh


// Unit headers
#include <utility/options/FileVectorOption.fwd.hh>

// Package headers
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>

// Project headers
#include <utility/file/FileName.hh>

// ObjexxFCL headers
#include <ObjexxFCL/char.functions.hh>


namespace utility {
namespace options {


/// @brief Program file vector option class
class FileVectorOption :
	public VectorOption_T_< FileVectorOptionKey, file::FileName >
{


private: // Types


	typedef  VectorOption_T_< FileVectorOptionKey, file::FileName >  Super;
	typedef  file::FileName  FileName;


public: // Types


	using Super::default_value;
	using Super::def;
	using Super::value;
	using Super::operator ();


public: // Creation


	/// @brief Default constructor
	inline
	FileVectorOption()
	{}


	/// @brief Key + description constructor
	inline
	FileVectorOption(
		FileVectorOptionKey const & key_a,
		std::string const & description_a
	) :
		Super( key_a, description_a )
	{}


	/// @brief Clone this
	inline
	FileVectorOption *
	clone() const override
	{
		return new FileVectorOption( *this );
	}


	/// @brief Destructor
	inline

	~FileVectorOption() {}


public: // Methods


	/// @brief Default value assignment
	inline
	FileVectorOption &
	default_value( std::string const & value_a )
	{
		Super::default_value( FileName( value_a ) );
		return *this;
	}


	/// @brief Default value assignment
	inline
	FileVectorOption &
	def( std::string const & value_a )
	{
		Super::default_value( FileName( value_a ) );
		return *this;
	}


	/// @brief Value assignment
	inline
	FileVectorOption &
	value( std::string const & value_a )
	{
		Super::value( FileName( value_a ) );
		return *this;
	}


	/// @brief Value assignment
	inline
	FileVectorOption &
	operator ()( std::string const & value_a )
	{
		Super::operator ()( FileName( value_a ) );
		return *this;
	}


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
		return "(F" + size_constraint_string() + ')';
	}


protected: // Methods


	/// @brief Value of a string
	inline
	Value
	value_of( std::string const & value_str ) const override
	{
		return FileName( value_str );
	}


}; // FileVectorOption


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_FileVectorOption_HH
