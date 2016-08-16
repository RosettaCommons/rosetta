// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/membrane/mp_dock_setup.cc
///
/// @brief  RosettaMP Membrane Protein-Protein Docking Protocol Setup Application
/// @details Setup poses for 2-body docking
///
/// @author  Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// @note   Last Updated: 5/18/15

// App headers
#include <devel/init.hh>

// Project headers
#include <protocols/docking/membrane/MPDockingSetupMover.hh>
#include <core/pose/Pose.fwd.hh>

// Package headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "apps.public.membrane.mp_dock_setup" );

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

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
