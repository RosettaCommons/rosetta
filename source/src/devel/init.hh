// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/devel/init.hh
/// @brief  Initialization function to ensure all load-time factory registration occurs
///         for classes that live in the devel library.  devel::init() calls
///         protocols::init::init(), which in turn calls core::init::init()
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_devel_init_hh
#define INCLUDED_devel_init_hh

// Utility headers
#include <utility/vector1.fwd.hh>

// C++ headers
#include <string>

namespace devel {

/// @brief Command line init() version
void init( int argc, char * argv [] );

/// @brief Wrapper for the command line version
void init( utility::vector1< std::string > const & args );

}

#endif

