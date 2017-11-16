// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/membrane/mp_dock.cc
///
/// @brief  RosettaMP Membrane Protein-Protein Docking Protocol
/// @details Dock to wproteins in the membrane
///
/// @author  Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// @note   Last Updated: 5/18/15

// App headers
#include <devel/init.hh>

// Project headers
#include <protocols/docking/membrane/MPDockingMover.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "apps.public.membrane.mp_dock" );

int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols::jd2;
		using namespace protocols::docking::membrane;


		// initialize options, RNG, and factory-registrators
		devel::init(argc, argv);

		MPDockingMoverOP mpdm( new MPDockingMover() );
		JobDistributor::get_instance()->go(mpdm);
	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}

	return 0;
}
