// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    glycosyltransferase.cc
/// @brief   This application performs a simulated glycosylation of an input structure.
/// @author  Labonte <JWLabonte@jhu.edu>


// Project headers
#include <devel/init.hh>

#include <protocols/enzymatic_movers/GlycosyltransferaseMover.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>


using namespace std;
using namespace utility;
using namespace core;
using namespace protocols;


// Constants & Type Definitions ///////////////////////////////////////////////
int const SUCCESS( 0 );
int const FAILURE( -1 );


// Main ///////////////////////////////////////////////////////////////////////
int
main( int argc, char *argv[] )
{
	try {
		// Initialize Rosetta.
		cout << "Initializing Rosetta..." << endl;
		devel::init( argc, argv );

		// Construct the protocol.
		enzymatic_movers::GlycosyltransferaseMoverOP protocol( new enzymatic_movers::GlycosyltransferaseMover );

		// Distribute the mover.
		jd2::JobDistributor::get_instance()->go( protocol );
	} catch (utility::excn::Exception const & e ) {
		cerr << "Caught exception: " << e.msg() << endl;
		return FAILURE;
	}
	return SUCCESS;
}
