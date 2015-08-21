// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/init/init.hh
/// @brief  Initialization function to ensure all load-time factory registration occurs
///         for classes that live in the protocols library.  protocols::init::init() calls core::init::init(),
///         and devel::init() call protocols::init::init().
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_init_init_HH
#define INCLUDED_protocols_init_init_HH

// Utility headers
#include <utility/vector1.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace init {

/// @brief Command line init() version
void init( int argc, char * argv [] );

/// @brief Wrapper for the command line version
void init( utility::vector1< std::string > const & args );

} //namespace init
} //namespace protocols

#endif //INCLUDED_protocols_init_init_HH
