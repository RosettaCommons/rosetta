// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/AnyVectorOption.hh
/// @brief  Program any vector-valued option abstract base class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_AnyVectorOption_hh
#define INCLUDED_utility_options_AnyVectorOption_hh


// Unit headers
#include <utility/options/AnyVectorOption.fwd.hh>

// Package headers
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>


namespace utility {
namespace options {


/// @brief Program any vector-valued option abstract base class
template< typename T >
class AnyVectorOption :
	public VectorOption_T_< AnyVectorOptionKey, T >
{


private: // Types


	typedef  VectorOption_T_< AnyVectorOptionKey, T >  Super;


public: // Types


	// STL/boost style
	typedef  AnyVectorOptionKey  key_type;
	typedef  T  value_type;

	// Project style
	typedef  AnyVectorOptionKey  Key;
	typedef  T  Value;


protected: // Creation


	/// @brief Default constructor
	inline
	AnyVectorOption()
	= default;


	/// @brief Copy constructor
	inline
	AnyVectorOption( AnyVectorOption const & option ) :
		Super( option )
	{}


	/// @brief Key + description constructor
	inline
	AnyVectorOption(
		AnyVectorOptionKey const & key_a,
		std::string const & description_a
	) :
		Super( key_a, description_a )
	{}


public: // Creation


	/// @brief Clone this
	virtual
	AnyVectorOption *
	clone() const = 0;


	/// @brief Destructor
	inline
	virtual
	~AnyVectorOption() {}


protected: // Assignment


	/// @brief Copy assignment
	inline
	AnyVectorOption &
	operator =( AnyVectorOption const & option )
	{
		if ( this != &option ) {
			Super::operator =( option );
		}
		return *this;
	}


public: // Properties


	/// @brief Option type code string representation
	inline
	std::string
	type_string() const
	{
		return "(A" + Super::size_constraint_string() + ')';
	}


}; // AnyVectorOption


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_AnyVectorOption_HH
