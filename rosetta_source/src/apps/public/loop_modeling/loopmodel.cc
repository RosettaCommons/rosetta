// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief demo program for implementing loop relax + FA relax
/// @author Srivatsan Raman
/// @author James Thompson
/// @author Mike Tyka

#include <protocols/abinitio/AbrelaxApplication.hh>

#include <protocols/loops/LoopBuild.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <devel/init.hh>

#include <utility/vector1.hh>



////////////////////////////////////////////////////////
void *
LoopRelax_main_local( void* ) {
	 protocols::loops::LoopRelax_main( false );
	return 0;
}

////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	// options, random initialization
	protocols::abinitio::ClassicAbinitio::register_options();
	protocols::abinitio::AbrelaxApplication::register_options();
	devel::init( argc, argv );
	protocols::viewer::viewer_main( LoopRelax_main_local );
	return 0;
}

