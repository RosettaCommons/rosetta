// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/options/ScalarOption_T_.hh
/// @brief  Program scalar-valued option abstract base class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Modified by Sergey Lyskov

#ifdef USEMPI
#include <mpi.h> // Must go first
#endif
// Unit headers
#include <utility/options/mpi_stderr.hh>

// C++ headers
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <set>
#include <sstream>

namespace utility {
namespace options {


void mpi_safe_std_err( std::string msg ) {
#ifdef  __native_client__
	throw( std::string( msg ) ); // exceptions provide a good mechanism to get the rror message back to where the browser can display it. stdout/stderr are useless here.
#endif

	int mpi_rank( 0 );
#ifdef USEMPI
	/// Give a different RNG to each processor

	MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
#endif
	if ( mpi_rank == 0 ) {
		std::cerr << msg << std::endl;std::cerr.flush();
	}
}

}
}

