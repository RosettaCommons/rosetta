// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/revold.cc
/// @brief  An mpi compatible evolutionary algorithm for ligand design
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

// project headers
#include <protocols/ligand_evolution/EvolutionManager.hh>

// utility headers
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>

// C/C++ headers
#ifdef USEMPI
#include <mpi.h>
#endif

static basic::Tracer TR( "apps.revold" ); // NOLINT(cert-err58-cpp)

using namespace protocols::ligand_evolution;


int main( int argc, char* argv[] ) {

	try {

		devel::init(argc, argv);

		int rank = 0;
		int size = 1;
#ifdef USEMPI
        MPI_Init (&argc, &argv);
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
        MPI_Comm_size( MPI_COMM_WORLD, &size );

        if ( size < 10 ) {
            TR.Warning << "Running REvoLd less than 10 CPUs causes very large runtimes. Consider using more CPUs." << std::endl;
        }
#else
		TR.Warning << "Running REvoLd without MPI causes very large runtimes. Consider using MPI." << std::endl;
#endif

		EvolutionManager manager( rank );

		manager.init();
		manager.run( size );


#ifdef USEMPI
        MPI_Barrier( MPI_COMM_WORLD );
        MPI_Finalize();
#endif

		return 0;
	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
}

