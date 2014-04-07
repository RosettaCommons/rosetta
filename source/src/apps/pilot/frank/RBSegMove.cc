// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Srivatsan Raman
/// @author Frank DiMaio
#include <protocols/rbsegment_relaxRelax_main.hh>
#include <protocols/viewer/viewers.hh>
#include <devel/init.hh>

#include <basic/options/option_macros.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
// Auto-header: duplicate removed #include <devel/init.hh>
#include <basic/options/option.hh>

#include <utility/excn/Exceptions.hh>


// options
OPT_1GRP_KEY( Boolean, fpd, viewer )

////////////////////////////////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////

void *
RBSegmentRelax_local_main( void* ) {
	protocols::RBSegmentRelax_main( false );
	return 0;
}

int
main( int argc, char * argv [] )
{
    try {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	NEW_OPT( fpd::viewer, "viewer?", false );

	// options, random initialization
	devel::init( argc, argv );

	if ( option[ fpd::viewer ]() )
		protocols::viewer::viewer_main( RBSegmentRelax_local_main );
	else
		RBSegmentRelax_local_main ((void*)0);
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
    }
