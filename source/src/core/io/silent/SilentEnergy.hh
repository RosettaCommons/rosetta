// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/SilentEnergy.hh
///
/// @brief simple class for managing energies in SilentStruct objects.
/// @author James Thompson

#ifndef INCLUDED_core_io_silent_SilentEnergy_hh
#define INCLUDED_core_io_silent_SilentEnergy_hh

#include <core/types.hh>
#include <string>

namespace core {
namespace io {
namespace silent {

/// @brief Helper class for silent-file classes to keep track of energy information.
class SilentEnergy {
public:
	SilentEnergy():
		name_  (  "" ),
		value_ ( 0.0 ),
		weight_( 0.0 ),
		width_ (   8 ),
		string_value_( "" )
	{}

	SilentEnergy(
		std::string name,
		core::Real  value,
		core::Real  weight,
		int         width
	) :
		name_  ( name  ),
		value_ ( value ),
		weight_( weight),
		width_ ( width ),
		string_value_( "" )
	{}

	SilentEnergy(
		std::string name,
		std::string value,
		int         width
	) :
		name_  ( name  ),
		value_ ( 0 ),
		weight_( 1.0 ),
		width_ ( width ),
		string_value_( value )
	{}

	std::string const& name() const {
		return name_;
	}

	void name( std::string const & new_name ) {
		name_ = new_name;
	}

	core::Real value() const {
		return value_;
	}

	void value( core::Real const & new_value ) {
		value_ = new_value;
	}

	std::string const& string_value() const {
		return string_value_;
	}

	void value( std::string const & new_string_value ) {
		string_value_ = new_string_value;
	}

	core::Real weight() const {
		return weight_;
	}

	void weight( core::Real const & new_weight ) {
		weight_ = new_weight;
	}

	int width() const {
		return width_;
	}

	void width( int const & new_width ) {
		width_ = new_width;
	}

private:
	std::string name_;

	core::Real value_;
	core::Real weight_;

	int width_;
	std::string string_value_;
}; // class SilentEnergy

} // namespace silent
} // namespace io
} // namespace core

#endif
