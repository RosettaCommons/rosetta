// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/VectorOption.hh
/// @brief  Program vector-valued option interface class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_options_VectorOption_hh
#define INCLUDED_utility_options_VectorOption_hh


// Unit headers
#include <utility/options/VectorOption.fwd.hh>

// Package headers
#include <utility/options/Option.hh>


namespace utility {
namespace options {


/// @brief Program vector-valued option interface class
class VectorOption :
	public Option
{


private: // Types


	typedef  Option  Super;


protected: // Creation


	/// @brief Default constructor
	inline
	VectorOption()
	{}


	/// @brief Copy constructor
	inline
	VectorOption( VectorOption const & ) = default;


public: // Creation


	/// @brief Clone this

	VectorOption *
	clone() const override = 0;


	/// @brief Destructor
	inline

	~VectorOption() override
	= default;


protected: // Assignment


	/// @brief Copy assignment
	inline
	VectorOption &
	operator =( VectorOption const & )
	{
		return *this;
	}


public: // Methods


	/// @brief Activate

	VectorOption &
	activate() override = 0;


	/// @brief Deactivate

	VectorOption &
	deactivate() override = 0;


	/// @brief Set to default value, if any

	VectorOption &
	to_default() override = 0;


	/// @brief Clear

	VectorOption &
	clear() override = 0;


	/// @brief Value assignment from a command line string

	VectorOption &
	cl_value( std::string const & value_str ) override = 0;

	/// @brief Value assignemt from a command line string but without
	/// a VectorOption & return type. This will separate arguments into
	/// blocks grouped by quotes, and then separate the non-quote-delimited
	/// arguments by whitespace

	void
	set_cl_value( std::string const & value_str ) override {
		cl_value( value_str );
	}

	/// @brief Fixed number of values required assignment
	virtual
	VectorOption &
	n( Size const n_a ) = 0;


	/// @brief Lower number of values allowed assignment
	virtual
	VectorOption &
	n_lower( Size const n_a ) = 0;


	/// @brief Upper number of values allowed assignment
	virtual
	VectorOption &
	n_upper( Size const n_a ) = 0;


public: // Properties


	/// @brief Fixed number of values required?
	virtual
	bool
	fixed_size() const = 0;


	/// @brief Fixed number of values required (zero if none)
	virtual
	Size
	n() const = 0;


	/// @brief Lower number of values allowed (zero if none)
	virtual
	Size
	n_lower() const = 0;


	/// @brief Upper number of values allowed (zero if none)
	virtual
	Size
	n_upper() const = 0;


	/// @brief Legal or inactive default value?
	virtual
	bool
	legal_default_value() const = 0;


	/// @brief Legal default value size?
	virtual
	bool
	legal_default_size() const = 0;


	/// @brief Legal value?
	virtual
	bool
	legal_value() const = 0;


	/// @brief Legal value size?
	virtual
	bool
	legal_size() const = 0;


}; // VectorOption


} // namespace options
} // namespace utility


#endif // INCLUDED_utility_options_VectorOption_HH
