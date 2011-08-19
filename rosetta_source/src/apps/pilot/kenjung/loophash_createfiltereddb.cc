// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Mike Tyka
/// @author Ken Jung
/// @brief

// libRosetta headers
//#include <basic/options/keys/in.OptionKeys.gen.hh>
//j#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <basic/Tracer.hh>
// C++ headers
#include <iostream>
#include <utility/excn/Exceptions.hh>

#include <core/init.hh>

#ifdef USEMPI
#include <mpi.h>
#endif

using basic::T;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;
using namespace protocols::loophash;
static basic::Tracer TR("main");


// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif


int
main( int argc, char * argv [] )
{


	// initialize core
	core::init(argc, argv);

	// Shouldn't these ints be core::size?
	// Default num_partitions = 1
	int num_partitions = option[lh::num_partitions]();
	int assigned_num = 1;
// No point in doing this without MPI
#ifndef USEMPI
	std::cerr << "You cannot use loophash_createfiltereddb without MPI!" << std::endl; 
	std::cout << "You cannot use loophash_createfiltereddb without MPI!" << std::endl; 
#endif

#ifdef USEMPI
	int mpi_rank_, mpi_npes_;
	MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &mpi_rank_ ) );
	MPI_Comm_size( MPI_COMM_WORLD, ( int* )( &mpi_npes_ ) );

	// unless you are rank under number of partitions, do nothing
	if( mpi_rank_ < num_partitions ) {
        assigned_num = mpi_rank_;
        
        utility::vector1 < core::Size > loop_sizes = option[lh::loopsizes]();
        LoopHashLibraryOP loop_hash_library = new LoopHashLibrary( loop_sizes, num_partitions, assigned_num );
        
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

        // Now master concats them all, throwing out nonunique fragments
        TR << "Finished my partition of the loophash library, waiting until all partitions are made" << std::endl;
        MPI_Barrier( MPI_COMM_WORLD );
        TR << "Restarting work since everyone's done with their work" << std::endl;
        if ( mpi_rank_ != 0 ) {
            // Delete libraries to regain heap memory
            // Since they're done anyways, do I need to do this?
            //loop_hash_library.reset_to_null();  // I think this is causing a segfault for some reason
            TR << "Not master, so no more work for me." << std::endl;
        } else {
            // delete old library and create a new empty library as the master library
            LoopHashLibraryOP loop_hash_library = new LoopHashLibrary( loop_sizes, 1, 0);
            utility::vector1< core::Real > rms_cutoff = option[lh::createdb_rms_cutoff]();
            if( rms_cutoff.size() != loop_sizes.size() ) {
                TR << "Number of RMS cutoff values does not match number of loop sizes" << std::endl;
                MPI_Abort( MPI_COMM_WORLD, 1 ); 
            }
            LoopHashLibraryOP second_loopdb;
            for( int i = 0; i < num_partitions ; i++) {
                second_loopdb = new LoopHashLibrary( loop_sizes, num_partitions, i );
                // Load that partition
                try {
                    second_loopdb->load_db();
										// Merge with master db throwing out uniques
										TR << "Merging master db with partition db " << i+1 << std::endl;
										loop_hash_library->merge( second_loopdb, rms_cutoff );
										TR << "Merge with partition " << i+1 << " finished." << std::endl;
										// Delete intermediate partition file
										second_loopdb->delete_db();
										// Future: Extra step that strips out backbone sequences that aren't used
										// Would be a pain because all the indexes would be shifted
                } catch ( utility::excn::EXCN_Base& excn ) {
                    excn.show( TR );
                    MPI_Abort( MPI_COMM_WORLD, 1 );
                }
            }
            TR << "All partitions merged, sorting and saving partition" << std::endl;
            //Save master db ( can add multiple partition saving in future )
            loop_hash_library->sort();
						try {
								loop_hash_library->save_db();
						} catch ( utility::excn::EXCN_Base& excn ) {
								excn.show( TR );
								MPI_Abort( MPI_COMM_WORLD, 1 );
						}
            TR << "Filtered DB successfully created." << std::endl;
       }
   }
   MPI_Finalize();
#endif
	return 0;
}
