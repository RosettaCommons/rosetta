// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers

#include <protocols/loops/Loops.hh>
#include <devel/init.hh>

// option key includes


//Auto Headers
#include <basic/options/keys/OptionKeys.hh>
#include <utility/excn/Exceptions.hh>


int
main( int argc, char * argv [] )
{
    try {
    	devel::init(argc, argv);
    	using namespace basic::options;
    	using namespace basic::options::OptionKeys;

    	// parse loops file Loops
    	//protocols::loops::Loops loops;
    	//std::string filename( option[ OptionKeys::loops::loop_file ]().name() );
    	//loops.read_loop_file( filename );   // <== TODO: select these using density score
    		protocols::loops::Loops myloops;
    		myloops.read_loop_file("testloops.loopfile" );
    		myloops.grow_all_loops(202, 4.0 );
    } catch (utility::excn::Exception const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}

