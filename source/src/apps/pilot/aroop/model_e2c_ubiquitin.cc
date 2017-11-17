// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Aroop Sircar ( aroopsircar@yahoo.com )
/// @brief

#include <protocols/jobdist/standard_mains.hh>

// Rosetta Headers
#include <basic/Tracer.hh>
#include <devel/init.hh>

#include <protocols/antibody_legacy/Ubiquitin_E2C_Modeler.hh>
#include <protocols/moves/Mover.hh>

////////////////////////////////////////////////////////
using basic::Error;
using basic::Warning;

using namespace core;

int
main( int argc, char * argv [] )
{
	try {

		using namespace protocols;
		using namespace protocols::jobdist;
		using namespace protocols::moves;

		// initialize core
		devel::init( argc, argv );

		MoverOP e2ubdock = new ub_e2c::ubi_e2c_modeler();
		protocols::jobdist::main_plain_mover( *e2ubdock );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}


