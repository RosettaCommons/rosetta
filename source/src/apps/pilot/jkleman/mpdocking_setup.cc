// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       mpdocking_setup.cc
/// @brief      Uses MPDockingSetup mover
///    CURRENTLY ONLY WORKS FOR 2 POSES!!!
/// @author     JKLeman (julia.koehler1982@gmail.com)
/// @note       Last Modified (10/16/14)

// App headers
#include <devel/init.hh>

// Project headers
#include <protocols/docking/membrane/MPDockingSetupMover.hh>
#include <core/pose/Pose.fwd.hh>

// Package headers
#include <basic/Tracer.hh>

static basic::Tracer TR( "apps.pilot.jkleman.mpdocking_setup" );

////////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {

		using namespace protocols::docking::membrane;

		// initialize options, RNG, and factory-registrators
		devel::init(argc, argv);

		// create empty pose
		core::pose::Pose pose;

		// create MPdockingSetup mover and apply to pose
		MPDockingSetupMoverOP mpdsm( new MPDockingSetupMover() );
		mpdsm->apply(pose);

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
