// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file   apps/pilot/ronj/optE_main.cc
/// @brief  Main function for optimizing weights for the surface term
/// @author Ron Jacak (ron.jacak@gmail.com)


/// MPI
#ifdef USEMPI
#include <mpi.h>
#endif

/// Core headers
#include <devel/init.hh>

// Protocol headers
#include <protocols/optimize_weights/IterativeOptEDriver.hh>

int
main( int argc, char * argv [] ) {

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


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
