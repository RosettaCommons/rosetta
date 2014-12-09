// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Mike Tyka
/// @brief

// libRosetta headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
// C++ headers
#include <iostream>
// option key includes

#include <basic/options/keys/lh.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


#ifdef USEMPI
#include <mpi.h>
#endif

static thread_local basic::Tracer TR( "main" );


// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif


int
main( int argc, char * argv [] )
{
    try {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;
	using namespace protocols;


	// initialize core
	devel::init(argc, argv);

   // Shouldn't these ints be core::size?
   // Default num_partitions = 1
    int num_partitions = option[lh::num_partitions]();
    int assigned_num = 1;

#ifdef USEMPI
	int mpi_rank_, mpi_npes_;
	MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &mpi_rank_ ) );
	MPI_Comm_size( MPI_COMM_WORLD, ( int* )( &mpi_npes_ ) );

	// unless you are rank one - go into infinite sleep loop
	if( mpi_rank_ > num_partitions - 1 ){
		TR << "NOT UNDER RANK " << num_partitions << ": Sleeping .. " << std::endl;
		while(true){
			sleep( 10 );
		}
	} else {
        assigned_num = mpi_rank_;
    }
#endif
	utility::vector1 < core::Size > loop_sizes = option[lh::loopsizes]();
	loophash::LoopHashLibraryOP loop_hash_library( new loophash::LoopHashLibrary( loop_sizes, num_partitions, assigned_num ) );
	loop_hash_library->create_db();
	loop_hash_library->save_db();
	TR << "Finished creating loophash library" << std::endl;
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}









