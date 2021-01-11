// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/jackmaguire/test_rosetta_thread_manager_preallocation.cc
/// @brief An application to test the RosettaThreadManager's preallocation API
/// @author Jack Maguire, jackmaguire1444@gmail.com

// devel headers
#include <devel/init.hh>

// core headers
#include <core/types.hh>

// utility headers
#include <utility/excn/Exceptions.hh>
#include <numeric/random/random.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/multithreading.OptionKeys.gen.hh>

// multithreading headers
#ifdef MULTI_THREADED
#include <basic/thread_manager/RosettaThreadManager.hh>
#include <basic/thread_manager/RosettaThreadPool.hh>
#include <basic/thread_manager/RosettaThreadAssignmentInfo.hh>
#endif //MULTI_THREADED

#include <chrono>

static basic::Tracer TR( "apps.pilot.jackmaguire.test_rosetta_thread_manager_preallocation" );

struct SleepForATenthOfASecond {
	void
	operator()(){
		std::chrono::duration< double > elapsed_seconds;
		auto start = std::chrono::steady_clock::now();
		do {
			auto end = std::chrono::steady_clock::now();
			elapsed_seconds = end - start;
		} while( elapsed_seconds.count() < 0.1 ); //This is busy waiting - not a recommended way to sleep
	}
};

//double
//round( double d ){
// return double(int( d * 10 )) / 10;
//}

int
main( int argc, char * argv [] )
{
	try {

		devel::init( argc, argv );

#ifdef MULTI_THREADED

		using namespace basic::thread_manager;

		RosettaThreadManager & manager = * RosettaThreadManager::get_instance();
		RosettaThreadAssignmentInfo assignment_info( RosettaThreadRequestOriginatingLevel::APPLICATIONS_OR_APPLICATION_PROTOCOLS );

		core::Size const nthreads = 3;
		RosettaThreadAllocation allocation = manager.reserve_threads( nthreads, assignment_info );
		core::Size const actual_nthreads = allocation.thread_ids.size() + 1;

		if ( actual_nthreads == 1 ) {
			utility_exit_with_message( "test_rosetta_thread_manager_preallocation is being run with only run thread, which is a poor way to test multithreading." );
		}

		TR << "Asked for " << nthreads << " threads. Was given " << actual_nthreads << " threads." << std::endl;

		utility::vector1< SleepForATenthOfASecond > const work( actual_nthreads );
		auto start = std::chrono::steady_clock::now();
		manager.do_work_vector_in_threads_no_locking( work, allocation, assignment_info );
		auto end = std::chrono::steady_clock::now();
		std::chrono::duration< double > elapsed_seconds = end - start;

		core::Real const expected_serial_runtime = 0.1 * actual_nthreads; //seconds
		core::Real const expected_parallel_runtime = 0.1; //seconds
		core::Real const runtime = elapsed_seconds.count();

		runtime_assert( runtime < expected_serial_runtime );//Already asserting that actual_nthreads > 1

		TR << std::fixed << std::setprecision(2) << std::endl;
		TR << "expected_serial_runtime: " << expected_serial_runtime << " seconds" << std::endl;
		TR << "expected_parallel_runtime: " << expected_parallel_runtime << " seconds" << std::endl;
		TR << "runtime: " << runtime << " seconds" << std::endl;

#endif

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
