// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test_RingConformationMover.cc
/// @brief   Pilot application source code for testing RingConformationMover.
/// @author  Labonte <JWLabonte@jhu.edu>

// Unit headers
#include <protocols/simple_moves/RingConformationMover.hh>

// Project headers
#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>


int main( int argc, char *argv[] )
{
	using namespace std;
	using namespace protocols::simple_moves;
	using namespace protocols::jd2;

	try {
		// Initialize core.
		devel::init( argc, argv );

		// Construct the mover.
		RingConformationMoverOP my_mover( new RingConformationMover() );

		// Distribute the mover.
		JobDistributor::get_instance()->go( my_mover );
	} catch (utility::excn::Exception const & e ) {
		cerr << "Caught exception: " << e.msg() << endl;
		return -1;
	}
	return 0;
}
