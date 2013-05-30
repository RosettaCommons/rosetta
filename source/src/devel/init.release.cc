// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/devel/init.release.cc
/// @brief  Initialization function to ensure all load-time factory registration occurs
///         for classes that live in the devel library.  devel::init() calls
///         protocols::init::init(), which in turn calls core::init::init().  Devel library does
///         not exist in release; so the release version of this file is nearly empty.
/// @author Steven Lewis (smlewi@gmail.com)

#include <devel/init.hh>
#include <protocols/init/init.hh>

namespace devel {

void init( int argc, char * argv [] )
{
	protocols::init::init( argc, argv );
}

void init( utility::vector1< std::string > const & args )
{
	protocols::init::init( args );
} // init

} // devel
