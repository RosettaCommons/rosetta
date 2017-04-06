// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/ralford/enumerate_rotamers.cc
/// @brief Create PDB files for all possible rotamers of D, E, N, and Q
/// @author Rebecca Alford (ralford3@jhu.edu)

// COOL! Now - whats left to do - generalize the script, maybe add the mover to master and create a pull request
// may also want to print out the rotamer probabiities

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/simple_moves/DumpSingleResidueRotamers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

// C++ Headers
#include <cstdlib>
#include <string>

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.ralford.enumerate_rotamers" );

int
main( int argc, char * argv [] )
{

	using namespace protocols::moves;
	using namespace protocols::simple_moves;

	try {
		devel::init( argc, argv );
		DumpSingleResidueRotamersOP dump_rotamers( new DumpSingleResidueRotamers() );
		protocols::jd2::JobDistributor::get_instance()->go( dump_rotamers );
		return 0;
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "Caught Exception " << e.msg() << std::endl;
		return -1;
	}

}
