// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/the_developer/my_app.cc
/// @brief Here is a breif description of my_app
/// @author the_developer (the_developer@foo.edu)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/jd2/JobDistributor.hh>

#include <protocols/moves/PyMOLMover.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>

static basic::Tracer TR( "apps.pilot.the_developer.my_app" );

int
main( int argc, char * argv [] )
{
	try {

		devel::init( argc, argv );

		protocols::moves::PyMOLMoverOP pmm( new protocols::moves::PyMOLMover() );
		pmm->keep_history( true );

		protocols::jd2::JobDistributor::get_instance()->go( pmm );

		TR << "************************************d**o**n**e***********************************" << std::endl;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
