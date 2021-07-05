// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/vmullig/test_rosetta_thread_manager_advanced_API.cc
/// @brief An application to test the Rosetta thread manager.  This application carries out three levels of
/// requests for threads.  Depending on the number of threads launched and the options provided, some might
/// not be possible to field, resulting in the thread manager assigning work to fewer threads than requested.
/// @details The highest-level threads work together to calcualte a times table.  Note that this tests the
/// advanced API for the RosettaThreadManager, which is not accessible to most Rosetta modules.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// devel headers
#include <devel/init.hh>

// protocol headers

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/multithreading.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>

// Multithreading includes

// STL includes

#include <core/types.hh> // AUTO IWYU For Size

#ifdef MULTI_THREADED
#include <basic/thread_manager/RosettaThreadManager.hh>
#include <basic/thread_manager/RosettaThreadAssignmentInfo.hh>
#include <basic/thread_manager/RosettaThreadManagerAdvancedAPIKey.hh>
#endif

#define NUMBER_OF_NUMBERS 1000
#define MAX_MULTIPLE 100

static basic::Tracer TR_main("test_rosetta_thread_manager_main");
static basic::Tracer TR_level0("test_rosetta_thread_manager_level_0");
static basic::Tracer TR_level1("test_rosetta_thread_manager_level_1");
static basic::Tracer TR_level2("test_rosetta_thread_manager_level_2");
static basic::Tracer TR_level3("test_rosetta_thread_manager_level_3");

OPT_KEY( Integer, level1_delay )
OPT_KEY( Integer, level2_delay )
OPT_KEY( Integer, level3_delay )
OPT_KEY( Integer, level1_threads )
OPT_KEY( Integer, level2_threads )
OPT_KEY( Integer, level3_threads )

/// @brief Storage for some data.
/// @note Since this is a dumb container class, all member variables are public.
class AppSettings {
public:

	core::Size nthreads_ = 0;

	core::Size delay_level1_ = 0;
	core::Size delay_level2_ = 0;
	core::Size delay_level3_ = 0;

	core::Size threads_level1_ = 0;
	core::Size threads_level2_ = 0;
	core::Size threads_level3_ = 0;

	/// @brief Initialization constructor.  Note that since this is a dumb container class,
	/// all member variables are public.
	AppSettings(
		core::Size const nthreads,
		core::Size const delay_level1,
		core::Size const delay_level2,
		core::Size const delay_level3,
		core::Size const threads_level1,
		core::Size const threads_level2,
		core::Size const threads_level3
	) :
		nthreads_( nthreads ),
		delay_level1_( delay_level1 ),
		delay_level2_( delay_level2 ),
		delay_level3_( delay_level3 ),
		threads_level1_( threads_level1 ),
		threads_level2_( threads_level2 ),
		threads_level3_( threads_level3 )
	{}
};

/// @brief Indicate required options.
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	NEW_OPT( level1_delay, "The lag time, in microseconds, for the task performed by the first-level threads before they launch threads.", 100 );
	NEW_OPT( level2_delay, "The lag time, in microseconds, for the task performed by the second-level threads before they launch threads.", 100 );
	NEW_OPT( level3_delay, "The lag time, in microseconds, for the task performed by the third-level threads before they launch threads.", 100 );
	NEW_OPT( level1_threads, "The number of level one worker threads that will be requested (including the master thread).", 4 );
	NEW_OPT( level2_threads, "The number of level two worker threads that will be requested from each level one thread (including the calling thread).", 4 );
	NEW_OPT( level3_threads, "The number of level three worker threads that will be requested from each level two thread (including the calling thread).", 4 );
	option.add_relevant( multithreading::total_threads );
}

#ifdef MULTI_THREADED

class RosettaThreadManagerTests {

public:

	RosettaThreadManagerTests() = default;
	RosettaThreadManagerTests( RosettaThreadManagerTests const & ) = default;
	~RosettaThreadManagerTests() = default;

	/// @brief The function executed by level three threads.
	/// @details Imagine that this were the main loop for a multithreaded version of a low-level Rosetta component, like the interaction
	/// graph calculator or the gradient vector calculator.
	void
	level_three_function(
		basic::thread_manager::RosettaThreadAssignmentInfo const & l3_thread_assignment_info,
		AppSettings const & appsettings,
		utility::vector1< core::Size > const & numbers,
		utility::vector1< std::mutex > & mutexes,
		utility::vector1< basic::thread_manager::AtomicBoolContainer > & completed,
		utility::vector1< utility::vector1< core::Size > > & multiples
	) {
		core::Size const thread_index( basic::thread_manager::RosettaThreadManager::get_instance()->get_rosetta_thread_index() );

		TR_level3 << "Level three thread " << thread_index << " reporting in.  " << appsettings.threads_level3_ << " threads were requested at this level.  " << l3_thread_assignment_info.get_assigned_total_thread_count() << " were assigned." << std::endl;
		std::this_thread::sleep_for( std::chrono::microseconds( appsettings.delay_level3_ ) );

		//Do something here: compute a times table.
		TR_level3 << "Preparing to operate on data in level three thread " << thread_index << "." << std::endl;
		for ( core::Size i(1), imax(numbers.size()); i<=imax; ++i ) { //Loop through all entries in vector.
			if ( completed[i].contained_bool_.load() ) continue;
			std::lock_guard< std::mutex > lock( mutexes[i] );
			if ( completed[i].contained_bool_.load() || !multiples[i].empty() ) continue; //This row of the times table has already been calculated.
			completed[i].contained_bool_ = true;
			TR_level3 << "Level three thread " << thread_index << " computing multiples of " << numbers[i] << "." << std::endl;
			multiples[i].resize( MAX_MULTIPLE );
			for ( core::Size j(1); j<=MAX_MULTIPLE; ++j ) {
				multiples[i][j] = numbers[i] * j;
			}
		}

		TR_level3 << "Level three thread " << thread_index << " signing off." << std::endl;
	}

	/// @brief The function executed by level two threads.
	/// @brief Imagine that this were the main loop for the threads of a multi-threaded mover.
	void
	level_two_function(
		basic::thread_manager::RosettaThreadAssignmentInfo const & l2_thread_assignment_info,
		AppSettings const & appsettings,
		utility::vector1< core::Size > const & numbers,
		utility::vector1< std::mutex > & mutexes,
		utility::vector1< basic::thread_manager::AtomicBoolContainer > & completed,
		utility::vector1< utility::vector1< core::Size > > & multiples
	) {
		TR_level2 << "Level two thread " << basic::thread_manager::RosettaThreadManager::get_instance()->get_rosetta_thread_index() << " reporting in.  " << appsettings.threads_level2_ << " threads were requested at this level.  " << l2_thread_assignment_info.get_assigned_total_thread_count() << " were assigned." << std::endl;
		std::this_thread::sleep_for( std::chrono::microseconds( appsettings.delay_level2_ ) );

		//Launch level three job:
		basic::thread_manager::RosettaThreadManager* thread_manager( basic::thread_manager::RosettaThreadManager::get_instance() );
		basic::thread_manager::RosettaThreadAssignmentInfo l3_thread_assignment_info( basic::thread_manager::RosettaThreadRequestOriginatingLevel::CORE_GENERIC ); //Pretend that this request comes from a core layer.
		basic::thread_manager::RosettaThreadFunction fxn( std::bind( &RosettaThreadManagerTests::level_three_function, this, std::cref(l3_thread_assignment_info), std::cref( appsettings ), std::cref( numbers ), std::ref( mutexes ), std::ref( completed ), std::ref( multiples ) ) );
		thread_manager->run_function_in_threads( fxn, appsettings.threads_level3_, basic::thread_manager::RosettaThreadManagerAdvancedAPIKey(), l3_thread_assignment_info );

		TR_level2 << "Level two thread " << basic::thread_manager::RosettaThreadManager::get_instance()->get_rosetta_thread_index() << " signing off." << std::endl;
	}

	/// @brief The function executed by level one threads.
	/// @details Imagine that this were the main loop for, say, the jobs of a multi-threaded job distributor.
	void
	level_one_function(
		basic::thread_manager::RosettaThreadAssignmentInfo const & l1_thread_assignment_info,
		AppSettings const & appsettings,
		utility::vector1< core::Size > const & numbers,
		utility::vector1< std::mutex > & mutexes,
		utility::vector1< basic::thread_manager::AtomicBoolContainer > & completed,
		utility::vector1< utility::vector1< core::Size > > & multiples
	) {
		TR_level1 << "Level one thread " << basic::thread_manager::RosettaThreadManager::get_instance()->get_rosetta_thread_index() << " reporting in.  " << appsettings.threads_level1_ << " threads were requested at this level.  " << l1_thread_assignment_info.get_assigned_total_thread_count() << " were assigned." << std::endl;
		std::this_thread::sleep_for( std::chrono::microseconds( appsettings.delay_level1_ ) );

		//Launch level two job:
		basic::thread_manager::RosettaThreadManager* thread_manager( basic::thread_manager::RosettaThreadManager::get_instance() );
		basic::thread_manager::RosettaThreadAssignmentInfo l2_thread_assignment_info( basic::thread_manager::RosettaThreadRequestOriginatingLevel::PROTOCOLS_GENERIC ); //Pretend that this request comes from a protocols layer.
		basic::thread_manager::RosettaThreadFunction fxn( std::bind( &RosettaThreadManagerTests::level_two_function, this, std::cref(l2_thread_assignment_info), std::cref( appsettings ), std::cref( numbers ), std::ref( mutexes ), std::ref( completed ), std::ref( multiples )) );
		thread_manager->run_function_in_threads( fxn, appsettings.threads_level2_, basic::thread_manager::RosettaThreadManagerAdvancedAPIKey(), l2_thread_assignment_info );

		TR_level1 << "Level one thread " << basic::thread_manager::RosettaThreadManager::get_instance()->get_rosetta_thread_index() << " signing off." << std::endl;
	}

	/// @brief The function executed by the main thread.
	/// @details Imagine that this were the main function of a multi-threaded job distributor.
	void
	level_zero_function(
		AppSettings const & appsettings,
		utility::vector1< core::Size > const & numbers,
		utility::vector1< std::mutex > & mutexes,
		utility::vector1< basic::thread_manager::AtomicBoolContainer > & completed,
		utility::vector1< utility::vector1< core::Size > > & multiples
	) {
		// Create the thread manager and launch threads:
		TR_level0 << "Spinning up global Rosetta thread manager with " << appsettings.nthreads_ << " total threads." << std::endl;
		basic::thread_manager::RosettaThreadManager* thread_manager( basic::thread_manager::RosettaThreadManager::get_instance() );
		basic::thread_manager::RosettaThreadAssignmentInfo thread_assignment_info( basic::thread_manager::RosettaThreadRequestOriginatingLevel::APPLICATIONS_OR_APPLICATION_PROTOCOLS ); //This request comes from the applications layer.

		// Create the function with std::bind:
		basic::thread_manager::RosettaThreadFunction fxn( std::bind( &RosettaThreadManagerTests::level_one_function, this, std::cref(thread_assignment_info), std::cref( appsettings ), std::cref( numbers ), std::ref( mutexes ), std::ref( completed ),  std::ref( multiples ) ) );

		// Pass it to the Rosetta thread manager to run in threads.  The minimum will be 1 thread, and the maximum will be min(nthread_level1, nthread_total) threads.
		thread_manager->run_function_in_threads( fxn, appsettings.threads_level1_, basic::thread_manager::RosettaThreadManagerAdvancedAPIKey(), thread_assignment_info );

		TR_level0 << "All level 1 threads have completed.  Returning from level 0." << std::endl;
	}

};

#endif //MULTI_THREADED

/// @brief Entry point for program execution.
int
main( int argc, char * argv [] )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	try {

		register_options();
		devel::init( argc, argv );

		TR_main << "Starting test_rosetta_thread_manager_advanced_API." << std::endl;
		TR_main << "Pilot app created 14 May 2019 by Vikram K. Mulligan, Flatiron Institute." << std::endl;
		TR_main << "For questions, please e-mail vmulligan@flatironinstitute.org." << std::endl;

#ifdef MULTI_THREADED
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		// Fill an array with the numbers 1 to 1000.
		utility::vector1< core::Size > numbers( NUMBER_OF_NUMBERS );
		for ( core::Size i(1); i<=NUMBER_OF_NUMBERS; ++i ) { numbers[i] = i; }
		utility::vector1< std::mutex > mutexes( NUMBER_OF_NUMBERS );
		utility::vector1< basic::thread_manager::AtomicBoolContainer > completed( NUMBER_OF_NUMBERS ); //Auto-initializes to false.
		utility::vector1< utility::vector1< core::Size > > multiples( NUMBER_OF_NUMBERS );

		AppSettings appsettings (
			option[ multithreading::total_threads ](),
			option[ level1_delay ](),
			option[ level2_delay ](),
			option[ level3_delay ](),
			option[ level1_threads ](),
			option[ level2_threads ](),
			option[ level3_threads ]()
		);

		RosettaThreadManagerTests tests;

		tests.level_zero_function( appsettings, numbers, mutexes, completed,  multiples );

		TR_main << "Results:";
		for ( core::Size i(1); i<=NUMBER_OF_NUMBERS; ++i ) {
			TR_main << "\n";
			for ( core::Size j(1); j<=MAX_MULTIPLE; ++j ) {
				TR_main << multiples[i][j] << " ";
				runtime_assert( multiples[i][j] == i*j );
			}
		}
		TR_main << std::endl;

#else //not MULTI_THREADED

		utility_exit_with_message( "This application cannot be run in the single-threaded build of Rosetta.  Please use the \"-extras=cxx11thread\" build option when compiling." );

#endif //MULTI_THREADED

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	TR_main << "Finished test_rosetta_thread_manager_advanced_API.  Exiting with status 0 (no errors)." << std::endl;

	return 0;
}
