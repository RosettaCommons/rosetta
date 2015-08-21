// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/options/IntegerVectorOption.hh
/// @brief  Program integer vector option class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Modified by Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_utility_options_IntegerVectorOption_HH
#define INCLUDED_utility_options_IntegerVectorOption_HH


// Unit headers
#include <utility/options/IntegerVectorOption.fwd.hh>

// Package headers
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <cstdlib>
#include <iostream>


namespace utility {
namespace options {


/// @brief Program integer vector option class
class IntegerVectorOption :
	public VectorOption_T_< IntegerVectorOptionKey, int >
{


private: // Types


	typedef  VectorOption_T_< IntegerVectorOptionKey, int >  Super;


public: // Creation


	/// @brief Default constructor
	inline
	IntegerVectorOption()
	{}


	/// @brief Key + description constructor
	inline
	IntegerVectorOption(
		IntegerVectorOptionKey const & key_a,
		std::string const & description_a
	) :
		Super( key_a, description_a )
	{}


	/// @brief Clone this
	inline
	IntegerVectorOption *
	clone() const
	{
		return new IntegerVectorOption( *this );
	}


	/// @brief Destructor
	inline
	virtual
	~IntegerVectorOption()
	{}


public: // Properties


	/// @brief Is a string readable as this option's value type?
	inline
	bool
	is_value( std::string const & value_str ) const
	{
		return ObjexxFCL::is_int( value_str );
	}


	/// @brief Is a string readable as this option's value type and a legal command line value?
	inline
	bool
	is_cl_value( std::string const & value_str ) const
	{
		return is_value( value_str ) || ObjexxFCL::is_ints( value_str );
	}


	/// @brief Option type code string representation
	inline
	std::string
	type_string() const
	{
		return "(I" + size_constraint_string() + ')';
	}


protected: // Methods


	/// @brief Value of a string
	inline
	Value
	value_of( std::string const & value_str ) const
	{
		if ( ( value_str.empty() ) || ( ! ObjexxFCL::is_int( value_str ) ) ) {
			std::cerr << "ERROR: Illegal value for integer option -" << id()
				<< " specified: " << value_str << std::endl;
			std::exit( EXIT_FAILURE );
		}
		return ObjexxFCL::int_of( value_str );
	}

	/// @brief Value of a string
	inline
	Values
	values_of( std::string const & value_str ) const
	{
		if ( value_str.empty() || ( ! ObjexxFCL::is_ints( value_str ) ) ) {
			std::cerr << "ERROR: Illegal value for integer option -" << id()
				<< " specified: " << value_str << std::endl;
			std::exit( EXIT_FAILURE );
		}
		std::vector< int > std_vector_ints =  ObjexxFCL::ints_of( value_str );

		//convert to utility vector1. This is a workaround to prevent ObjexxFCL from knowing about vector1.
		Values vector1_ints;
		for ( Size n = 0; n < std_vector_ints.size(); n++ ) vector1_ints.push_back(  std_vector_ints[ n ] );
		return vector1_ints;
	}


}; // IntegerVectorOption


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_IntegerVectorOption_HH
