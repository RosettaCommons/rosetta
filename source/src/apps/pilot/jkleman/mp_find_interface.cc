// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    mp_find_interface.cc
/// @brief   Sample protein-protein interface in the membrane
/// @details last Modified: 4/4/15
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project headers
#include <protocols/docking/membrane/MPFindInterfaceMover.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "apps.pilot.jkleman.mp_find_interface" );

//////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols::jd2;
		using namespace protocols::docking::membrane;


		// initialize options, RNG, and factory-registrators
		devel::init(argc, argv);

		MPFindInterfaceMoverOP mfi( new MPFindInterfaceMover() );
		JobDistributor::get_instance()->go(mfi);
	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}

	return 0;
}
