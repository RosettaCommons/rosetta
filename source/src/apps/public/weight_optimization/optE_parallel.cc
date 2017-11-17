// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

/// MPI headers
#ifdef USEMPI
#include <mpi.h>
#endif

/// Project headers
#include <protocols/optimize_weights/IterativeOptEDriver.hh>
#include <devel/init.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>


int
main( int argc, char * argv [] )
{
	try {
#ifdef USEMPI
		MPI_Init(&argc, &argv);
#endif


		devel::init( argc, argv );

		protocols::optimize_weights::IterativeOptEDriver driver;
		driver.go();

#ifdef USEMPI
		MPI_Finalize();
#endif
	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
