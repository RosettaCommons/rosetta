// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/VectorOption.cc
/// @brief  Program vector-valued option interface class
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Unit headers
#include <utility/options/VectorOption.hh>

// Package headers
#include <utility/string_util.hh>

namespace utility {
namespace options {

/// @brief Value assignemt from a command line string but without
/// a VectorOption & return type.
/// This will separate arguments by whitespace and set the vector option accordingly.
/// If reset is true, the current values will be cleared first.
void
VectorOption::set_value( std::string const & value_str, bool reset ) {
	if ( reset ) {
		clear();
	}
	utility::vector1< std::string > ws_sep_values = utility::quoted_split( value_str );
	for ( std::string const & str : ws_sep_values ) {
		cl_value( str );
	}
}

} // namespace options
} // namespace utility

