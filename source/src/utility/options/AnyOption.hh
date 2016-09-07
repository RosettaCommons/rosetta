// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/AnyOption.hh
/// @brief  Program any scalar-valued option abstract base class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_AnyOption_hh
#define INCLUDED_utility_options_AnyOption_hh


// Unit headers
#include <utility/options/AnyOption.fwd.hh>

// Package headers
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/keys/AnyOptionKey.hh>


namespace utility {
namespace options {


/// @brief Program any scalar-valued option abstract base class
template< typename T >
class AnyOption :
	public ScalarOption_T_< AnyOptionKey, T >
{


private: // Types


	typedef  ScalarOption_T_< AnyOptionKey, T >  Super;


public: // Types


	// STL/boost style
	typedef  AnyOptionKey  key_type;
	typedef  T  value_type;

	// Project style
	typedef  AnyOptionKey  Key;
	typedef  T  Value;


protected: // Creation


	/// @brief Default constructor
	inline
	AnyOption()
	= default;


	/// @brief Copy constructor
	inline
	AnyOption( AnyOption const & option ) :
		Super( option )
	{}


	/// @brief Key + description constructor
	inline
	AnyOption(
		AnyOptionKey const & key_a,
		std::string const & description_a
	) :
		Super( key_a, description_a )
	{}


public: // Creation


	/// @brief Clone this
	virtual
	AnyOption *
	clone() const = 0;


	/// @brief Destructor
	inline
	virtual
	~AnyOption()
	= default;


protected: // Assignment


	/// @brief Copy assignment
	inline
	AnyOption &
	operator =( AnyOption const & option )
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
		return "A";
	}


}; // AnyOption


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_AnyOption_HH
