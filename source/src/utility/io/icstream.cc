// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/io/icstream.cc
/// @brief  Input channel stream wrapper base class,
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


// Unit header
#include <utility/io/icstream.hh>

// C++ headers
#include <iostream>


namespace utility {
namespace io {
namespace ic { // Predefined icstreams


/// @brief Wrapper around std::cin
icstream cin( std::cin );


} // namespace ic
} // namespace io
} // namespace utility
