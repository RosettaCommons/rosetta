// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief demo program for implementing loop relax + FA relax
/// @author Srivatsan Raman
/// @author James Thompson
/// @author Mike Tyka

#include <protocols/abinitio/AbrelaxApplication.hh>

#include <protocols/loop_build/LoopBuild.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <devel/init.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

static basic::Tracer TR( "apps.public.loops.loopmodel" );

////////////////////////////////////////////////////////
void *
LoopBuild_main_local( void* ) {
	protocols::loop_build::LoopBuild_main( false );
	return 0;
}

////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		// options, random initialization
		protocols::abinitio::ClassicAbinitio::register_options();
		protocols::abinitio::AbrelaxApplication::register_options();
		devel::init( argc, argv );

		// Old versions of loopmodel used -loops:input_pdb for PDB input.
		// Catch people using old versions of command lines and tell them to update.
		if ( basic::options::option[ basic::options::OptionKeys::loops::input_pdb ].user() ) {
			TR.Error << "loopmodel no longer uses the -loops:input_pdb flag -- consult recent documentation for the updated flags to use." << std::endl;
			utility_exit_with_message("Incorrect flag passed: -loops:input_pdb");
		}

		protocols::viewer::viewer_main( LoopBuild_main_local );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

