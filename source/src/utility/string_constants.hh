// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    utility/string_constants.hh
/// @brief   Commonly used string constants.
/// @author  Labonte <JWLabonte@jhu.edu>


// C++ Header
#include <string>

namespace utility {

static std::string const UPPERCASE_LETTERS( "ABCDEFGHIJKLMNOPQRSTUVWXYZ" );
static std::string const LOWERCASE_LETTERS( "abcdefghijklmnopqrstuvwxyz" );
static std::string const LETTERS( UPPERCASE_LETTERS + LOWERCASE_LETTERS );

static std::string const NUMERALS( "0123456789" );

static std::string const UPPERCASE_ALPHANUMERICS( UPPERCASE_LETTERS + NUMERALS );
static std::string const LOWERCASE_ALPHANUMERICS( LOWERCASE_LETTERS + NUMERALS );
static std::string const ALPHANUMERICS( LETTERS + NUMERALS );

}  // namespace utility
