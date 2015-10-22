// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/options/ResidueChainVectorOption.hh
/// @brief  Program integer vector option class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Modified by Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_utility_options_ResidueChainVectorOption_HH
#define INCLUDED_utility_options_ResidueChainVectorOption_HH

// Unit headers
#include <utility/options/ResidueChainVectorOption.fwd.hh>

// Package headers
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/keys/ResidueChainVectorOptionKey.hh>

// ObjexxFCL headers
#include <ObjexxFCL/char.functions.hh>

// C++ headers
//#include <cstdlib>
//#include <iostream>


// A mod of IntegerVectorOption that, by default, returns list of integers but more generally stores and
// returns chain information as well. Converts tags like A:1-4 B:1-3 into a pair of vectors,  [1 2 3 4 1 2 3], [A A A A B B B]


namespace utility {
namespace options {


/// @brief Program integer vector option class
class ResidueChainVectorOption :
	public VectorOption_T_< ResidueChainVectorOptionKey, int >
{


private: // Types

	typedef  VectorOption_T_< ResidueChainVectorOptionKey, int >  Super;

public: // Creation


	/// @brief Default constructor
	inline
	ResidueChainVectorOption()
	{}


	/// @brief Key + description constructor
	inline
	ResidueChainVectorOption(
		ResidueChainVectorOptionKey const & key_a,
		std::string const & description_a
	) :
		Super( key_a, description_a )
	{}


	/// @brief Clone this
	inline
	ResidueChainVectorOption *
	clone() const
	{
		return new ResidueChainVectorOption( *this );
	}


	/// @brief Destructor
	inline
	virtual
	~ResidueChainVectorOption()
	{}


public: // Properties


	/// @brief Is a string readable as this option's value type?
	inline
	bool
	is_value( std::string const & ) const {
		return true;
	}

	/// @brief Is a string readable as this option's value type and a legal command line value?
	inline
	bool
	is_cl_value( std::string const & value_str ) const {
		return ( ( value_str.empty() ) || ( ! ObjexxFCL::is_any_of( value_str[ 0 ], "-@" ) ) );
	}

	/// @brief Option type code string representation
	inline
	std::string
	type_string() const {
		return "(RC" + size_constraint_string() + ')';
	}

	// @brief Specialized function that converts tags like A:1-4 B:1-3 into a pair of  [1 2 3 4 1 2 3], [A A A A B B B]
	std::pair< utility::vector1<int>, utility::vector1<char> >
	resnum_and_chain() const;

	// @brief specialized cl_value operator that saves value_string.
	VectorOption_T_< ResidueChainVectorOptionKey, int > &
	cl_value( std::string const & value_str );

protected: // Methods

	/// @brief Value of a string
	Value
	value_of( std::string const & value_str ) const;

	/// @brief Value of a string
	Values
	values_of( std::string const & value_str ) const;


private:

	utility::vector1< std::string > value_strings_;

}; // ResidueChainVectorOption


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_ResidueChainVectorOption_HH
