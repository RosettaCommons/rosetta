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

#include <core/types.hh>
#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/membrane.OptionKeys.gen.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/abinitio/AbrelaxApplication.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

void*
my_main( void *)
{
	protocols::abinitio::AbrelaxApplication abrelax;
	abrelax.run();
	return 0 ;
}

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init( argc, argv );
		protocols::abinitio::AbrelaxApplication::register_options();
		protocols::abinitio::AbrelaxApplication abrelax;

		bool const view( option[ OptionKeys::membrane::view ] );

		if ( view ) {
			std::cout << "Start viewer mode " << std::endl;
			protocols::viewer::viewer_main( my_main );
		} else {
			abrelax.run();
		}
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

