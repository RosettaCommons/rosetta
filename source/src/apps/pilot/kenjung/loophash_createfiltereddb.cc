// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Mike Tyka
/// @author Ken Jung
/// @brief

// libRosetta headers
//#include <basic/options/keys/in.OptionKeys.gen.hh>
//j#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#ifdef USEMPI
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <utility/excn/Exceptions.hh>
#endif

#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

#include <devel/init.hh>

#include <core/types.hh>
#include <protocols/loophash/LoopHashSampler.fwd.hh>
#include <utility/vector1.hh>


#ifdef USEMPI
#include <mpi.h>
#endif

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;
using namespace protocols::loophash;
static basic::Tracer TR( "main" );


// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif


int
main( int argc, char * argv [] )
{
	try {

		// initialize core
		devel::init(argc, argv);

		// Shouldn't these ints be core::size?
		// Default num_partitions = 1
		// No point in doing this without MPI
#ifndef USEMPI
		std::cerr << "You cannot use loophash_createfiltereddb without MPI!" << std::endl;
		std::cout << "You cannot use loophash_createfiltereddb without MPI!" << std::endl;
#endif

#ifdef USEMPI
	int num_partitions = option[lh::num_partitions]();
	int assigned_num = 1;

	int mpi_rank_, mpi_npes_;
	MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &mpi_rank_ ) );
	MPI_Comm_size( MPI_COMM_WORLD, ( int* )( &mpi_npes_ ) );
	long starttime = time(NULL);

	// unless you are rank under number of partitions, do nothing
	if( mpi_rank_ < num_partitions ) {
        assigned_num = mpi_rank_;

        utility::vector1 < core::Size > loop_sizes = option[lh::loopsizes]();
				utility::vector1< core::Real > rms_cutoff = option[lh::createdb_rms_cutoff]();
				if( rms_cutoff.size() != loop_sizes.size() ) {
						TR << "Number of RMS cutoff values does not match number of loop sizes" << std::endl;
						MPI_Abort( MPI_COMM_WORLD, 1 );
				}

        LoopHashLibraryOP loop_hash_library( new LoopHashLibrary( loop_sizes, num_partitions, assigned_num ) );


        TR << "Creating partition..." << std::endl;
				loop_hash_library->create_db();
        TR << "Saving partition..." << std::endl;
				try {
            loop_hash_library->save_db();
        } catch ( utility::excn::EXCN_Base& excn ) {
            excn.show( TR );
        		TR << "Exception occured!" << std::endl;
            MPI_Abort( MPI_COMM_WORLD, 1 );
        }

				// move the master (mpi_rank = 0) db to a nonpartitioned db name, and perform initial filtering
				if( mpi_rank_ == 0) {
					LoopHashLibraryOP old_master = loop_hash_library;
					loop_hash_library = LoopHashLibraryOP( new LoopHashLibrary( loop_sizes, 1, 0) );
					loop_hash_library->merge( old_master, rms_cutoff );
					old_master->delete_db();
				}

        TR << "Finished my partition of the loophash library, waiting until all partitions are made" << std::endl;
        MPI_Barrier( MPI_COMM_WORLD );
        TR << "Restarting work since everyone's done with their work" << std::endl;

				// divide and conquer merge!
				core::Size divisor = 1;
				while( (float)num_partitions / divisor > 1 ) {
					divisor = divisor * 2;
					if ( mpi_rank_ == 0) TR << "Divisor = " << divisor << std::endl;
					if ( mpi_rank_ % divisor != 0 || mpi_rank_ + (int) divisor/2 >= (int) num_partitions ) {
						MPI_Barrier( MPI_COMM_WORLD );
					} else {
						LoopHashLibraryOP second_loopdb( new LoopHashLibrary( loop_sizes, num_partitions, mpi_rank_ + divisor/2 ) );
						try {
								second_loopdb->load_db();
								// Merge while throwing out uniques
								TR << "Merging parts " << mpi_rank_ +1 << " with " << mpi_rank_ + 1+ divisor/2 << std::endl;
								loop_hash_library->merge( second_loopdb, rms_cutoff );
								//save new tmp db if it will be loaded as a second_loopdb next round
								if( mpi_rank_ % (divisor * 2) != 0 || mpi_rank_ + (int) divisor >= (int) num_partitions ) {
									loop_hash_library->save_db();
									TR << "Saving " << mpi_rank_ << std::endl;
								}
								// Future: Extra step that strips out backbone sequences that aren't used
								// Would be a pain because all the indexes would be shifted
						} catch ( utility::excn::EXCN_Base& excn ) {
								excn.show( TR );
								MPI_Abort( MPI_COMM_WORLD, 1 );
						}
						MPI_Barrier( MPI_COMM_WORLD );
						if ( mpi_rank_ == 0 ) TR << "FINISHED ROUND " << divisor << std::endl;
					}
				}
				if( mpi_rank_ == 0 ) {
					TR << "All partitions merged, sorting and saving partition" << std::endl;
					//Save master db ( can add multiple partition saving in future )
					loop_hash_library->sort();
					try {
							loop_hash_library->save_db();
					} catch ( utility::excn::EXCN_Base& excn ) {
							excn.show( TR );
							MPI_Abort( MPI_COMM_WORLD, 1 );
					}

				} else {
					// delete intermediate partition
					loop_hash_library->delete_db();
				}

   }
   MPI_Finalize();
	long endtime = time(NULL);
	if( mpi_rank_ == 0 ) {
		TR << "Filtered DB successfully created." << std::endl;
		TR << "Total time: " << (endtime - starttime) / 60 << " min" << std::endl;
	}
#endif
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
