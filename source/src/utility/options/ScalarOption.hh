// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/ScalarOption.hh
/// @brief  Program scalar-valued option interface class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_ScalarOption_hh
#define INCLUDED_utility_options_ScalarOption_hh


// Unit headers
#include <utility/options/ScalarOption.fwd.hh>

// Package headers
#include <utility/options/Option.hh>


namespace utility {
namespace options {


/// @brief Program scalar-valued option interface class
class ScalarOption :
	public Option
{


private: // Types


	typedef  Option  Super;


protected: // Creation


	/// @brief Default constructor
	inline
	ScalarOption()
	{}


	/// @brief Copy constructor
	inline
	ScalarOption( ScalarOption const & ) = default;


public: // Creation


	/// @brief Clone this

	ScalarOption *
	clone() const override = 0;


	/// @brief Destructor
	inline

	~ScalarOption() {}


protected: // Assignment


	/// @brief Copy assignment
	inline
	ScalarOption &
	operator =( ScalarOption const & )
	{
		return *this;
	}


public: // Methods


	/// @brief Activate

	ScalarOption &
	activate() override = 0;


	/// @brief Deactivate

	ScalarOption &
	deactivate() override = 0;


	/// @brief Set to default value, if any

	ScalarOption &
	to_default() override = 0;


	/// @brief Clear

	ScalarOption &
	clear() override = 0;


	/// @brief Value assignment from a command line string

	ScalarOption &
	cl_value( std::string const & value_str ) override = 0;

	/// @brief Value assignemt from a command line string but without
	/// a ScalarOption & return type.

	void
	set_cl_value( std::string const & value_str ) override {
		cl_value( value_str );
	}

}; // ScalarOption


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_ScalarOption_HH
