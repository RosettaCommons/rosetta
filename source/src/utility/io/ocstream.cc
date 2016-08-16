// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/io/ocstream.cc
/// @brief  Output channel stream wrapper base class,
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


// Unit header
#include <utility/io/ocstream.hh>

// C++ headers
#include <iostream>


namespace utility {
namespace io {
namespace oc { // Predefined ocstreams


/// @brief Wrapper around std::cout
ocstream cout( std::cout );

/// @brief Wrapper around std::cerr
ocstream cerr( std::cerr );

/// @brief Wrapper around std::clog
ocstream clog( std::clog );


} // namespace oc
} // namespace io
} // namespace utility
