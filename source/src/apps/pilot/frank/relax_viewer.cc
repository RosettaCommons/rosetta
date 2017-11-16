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

#include <utility/excn/Exceptions.hh>

#include <protocols/jobdist/standard_mains.hh>
#include <basic/options/option.hh>

#include <devel/init.hh>

// C++ headers
#include <iostream>
#include <string>

#include <basic/Tracer.hh>

// option key includes
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/util.hh>
#include <protocols/relax/relax_main.hh>
#include <protocols/viewer/viewers.hh>
#include <basic/options/option_macros.hh>

using basic::Error;
using basic::Warning;


using namespace core;
using namespace protocols;

using utility::vector1;


///////////////////////////////////////////////////////////////////////////////
void *
relax_main_local( void* ) {
	relax::Relax_main( false );
	return 0;
}

////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols::jobdist;
		using namespace protocols::moves;
		using namespace scoring;
		using namespace basic::options;


		relax::ClassicRelax::register_options();
		register_options_universal_main();
		option.add_relevant( OptionKeys::in::file::fullatom );
		option.add_relevant( OptionKeys::relax::sequence );
		devel::init(argc, argv);

		// use viewer if flag given
		protocols::viewer::viewer_main( relax_main_local );
		relax_main_local(NULL);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}


