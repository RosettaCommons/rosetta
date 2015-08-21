// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    mp_quick_relax.cc
/// @brief   Do a quick relax run for a membrane protein
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project headers
#include <protocols/membrane/MPQuickRelaxMover.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "apps.pilot.jkleman.mp_quick_relax" );

//////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols::jd2;
		using namespace protocols::membrane;


		// initialize options, RNG, and factory-registrators
		devel::init(argc, argv);

		MPQuickRelaxMoverOP mqr( new MPQuickRelaxMover() );
		mqr->add_membrane_again( true );
		JobDistributor::get_instance()->go(mqr);
	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}

	return 0;
}
