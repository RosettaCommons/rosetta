// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadManager.cc
/// @brief A manager that maintains a threadpool and handles requests for threads for multithreaded
/// execution of functions.  This allows multithreading at many different levels in the Rosetta
/// library hierarchy, from job-level parallel execution down to parallel computation of a score,
/// gradient vector, or interaction graph.
/// @note In single-threaded builds, this object still exists.  It accepts vectors of work and executes
/// them directly, in this case.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#include <basic/thread_manager/RosettaThreadManager.hh>
#include <basic/thread_manager/RosettaThreadPool.hh>
#include <basic/thread_manager/RosettaThreadAssignmentInfo.hh>

// Utility headers
#include <utility/pointer/memory.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/multithreading.OptionKeys.gen.hh>

// C++ headers
#include <string>

// Construct tracer.
static basic::Tracer TR( "basic.thread_manager.RosettaThreadManager" );

namespace basic {
namespace thread_manager {

// Private methods ////////////////////////////////////////////////////////////

/// @brief Empty constructor.
RosettaThreadManager::RosettaThreadManager()
#ifdef MULTI_THREADED
:
	nthreads_( basic::options::option[ basic::options::OptionKeys::multithreading::total_threads ]() )
#endif
{
#ifdef MULTI_THREADED
	if ( nthreads_ == 0 ) {
		nthreads_ = std::thread::hardware_concurrency();
		if ( nthreads_ == 0 ) { //Could not determine number of threads (should be rare).
			nthreads_ = 1;
		}
	}
#endif
}

/// @brief Destructor.  Non-empty, since threads must be spun down.
RosettaThreadManager::~RosettaThreadManager() {
#ifdef MULTI_THREADED
	if ( thread_pool_ != nullptr ) {
		//TR << "Forcing termination of threads that are still running." << std::endl;
		thread_pool_->force_terminate_threads();
	}
#endif
}

#ifdef MULTI_THREADED
/// @brief Creates the thread pool if it has not yet been created.  Safe to call repeatedly.  The
/// thread_pool_mutex_ must be locked before calling this function!
/// @details Accesses the global options system to determine the number of threads to launch.
void
RosettaThreadManager::create_thread_pool() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	RosettaThreadManagerInitializationTracker::get_instance()->mark_thread_manager_initialization_as_begun();

	//#ifndef NDEBUG
	bool const is_unlocked( thread_pool_mutex_.try_lock() );
	if ( is_unlocked ) {
		thread_pool_mutex_.unlock();
		utility_exit_with_message( "Error in ThreadPoolManager::create_thread_pool():  This function was called, but the thread_pool_mutex_ was not locked!" );
	}
	//#endif

	if ( thread_pool_ == nullptr ) {
		TR << "Creating a thread pool of " << nthreads_ << " threads." << std::endl;
		thread_pool_ = utility::pointer::make_shared< RosettaThreadPool >( nthreads_, RosettaThreadPoolInstantiationKey() );
		thread_id_to_rosetta_thread_index_ = thread_pool_->get_thread_id_to_rosetta_thread_index_map();
	} else {
		TR.Warning << "WARNING WARNING WARNING!  RosettaThreadManager::create_thread_pool() was called, but the thread pool has already been created!" << std::endl;
	}

	RosettaThreadManagerInitializationTracker::get_instance()->mark_thread_manager_as_initialized();
}

/// @brief Trigger launch of threads.
/// @details Does nothing if threads already launched.
void
RosettaThreadManager::launch_threads() {
	std::lock_guard< std::mutex > lock( thread_pool_mutex_ );
	if ( thread_pool_ == nullptr ) {
		create_thread_pool();
	}
}
#endif

// Public methods /////////////////////////////////////////////////////////////

/// @brief BASIC API THAT SHOULD BE USED IN MOST CIRCUMSTANCES.  Given a vector of functions that were bundled with
/// their arguments with std::bind, each of which can be executed in any order and each of which is safe to execute
/// in parallel with any other, run all of these in threads.
/// @details The bundled functions should be atomistic pieces of work.  They should be bundled with their arguments with
/// std::bind, and the arguments should include the place to store output (i.e. they should return void).  These functions
/// should not handle any mutexes themselves, but should ensure that they are operating only on memory locations that no
/// other functions in the vector are operating on.
void
RosettaThreadManager::do_work_vector_in_threads(
	utility::vector1< RosettaThreadFunctionOP > const & vector_of_work,
	platform::Size const requested_thread_count,
#ifdef MULTI_THREADED
	RosettaThreadAssignmentInfoOP thread_assignment /*= nullptr*/
#else
	RosettaThreadAssignmentInfoOP /*= nullptr*/
#endif
) {
	if ( vector_of_work.size() == 0 ) {
		TR.Warning << "A work vector of size zero was passed to the RosettaThreadManager!  Duly returning without doing anything." << std::endl;
		return;
	}
#ifdef MULTI_THREADED
	platform::Size const nonzero_requested_thread_count( requested_thread_count == 0 ? nthreads_ : requested_thread_count );
	platform::Size const actual_threads_to_request( std::min( vector_of_work.size(), nonzero_requested_thread_count ) );
	utility::vector1< utility::thread::ReadWriteMutex > mutexes( vector_of_work.size() );
	utility::vector1< bool > jobs_completed( vector_of_work.size(), false );

	RosettaThreadFunctionOP fxn( utility::pointer::make_shared< RosettaThreadFunction > ( std::bind( &RosettaThreadManager::work_vector_thread_function, this, std::cref( vector_of_work ), std::ref( mutexes ), std::ref( jobs_completed ) ) ) );
	run_function_in_threads( fxn, actual_threads_to_request, RosettaThreadManagerAdvancedAPIKey(), thread_assignment );
#else
	//In non-threaded builds, just iterate through and do all of the work.
	runtime_assert_string_msg( requested_thread_count == 0 || requested_thread_count == 1, "Error in RosettaThreadManager::do_work_vector_in_threads(): In non-threaded builds, only one thread may be requested.  Compile Rosetta with the \"extras=cxx11thread\" option to enable multi-threading." );
	for ( platform::Size i(1), imax(vector_of_work.size()); i<=imax; ++i ) {
		(*vector_of_work[i])(); //Do the ith piece of work.
	}
#endif
}

/// @brief VARIANT BASIC API THAT SHOULD BE USED WHERE THE BASIC API CAN'T BE USED.  Given a vector of vectors of functions
/// that were bundled with their arguments with std::bind, run all of these in threads.  In this case, the bundled functions
/// are in groups, where the individual functions within a group can run in any order (and are safe to run concurrently),
/// but the groups must be run sequentially.  This is useful when, for example, you have a bunch of calculations to do, and
/// then some finalization tasks to do after the calculations are done, and you don't want to re-request threads.
/// @details The bundled functions should be atomistic pieces of work.  They should be bundled with their arguments with
/// std::bind, and the arguments should include the place to store output (i.e. they should return void).  These functions
/// should not handle any mutexes themselves, but should ensure that they are operating only on memory locations that no
/// other functions in the vector are operating on.
/// @note Under the hood, this sets up appropriate mutexes and then calls run_function_in_threads() to do the work.  The work
/// is done concurrently in 1 <= actual count <= min( requested thread count, total thread count ) threads.  The function blocks
/// until all threads have finished their work, which means that the individual work units should be small, that the longest-running
/// work unit should be short compared to the total runtime, and that the number of work units should be much greater than the
/// number of threads requested.
void
RosettaThreadManager::do_multistage_work_vector_in_threads(
	utility::vector1< utility::vector1< RosettaThreadFunctionOP > > const & multistage_vector_of_work,
	platform::Size const requested_thread_count,
#ifdef MULTI_THREADED
	RosettaThreadAssignmentInfoOP thread_assignment /*= nullptr*/
#else
	RosettaThreadAssignmentInfoOP /*= nullptr*/
#endif
) {
	if ( multistage_vector_of_work.size() == 0 ) {
		TR.Warning << "A work vector of size zero was passed to the RosettaThreadManager!  Duly returning without doing anything." << std::endl;
		return;
	}
#ifdef MULTI_THREADED
	platform::Size max_batch_size(0);
	for ( platform::Size i(1), imax(multistage_vector_of_work.size()); i<=imax; ++i ) {
		platform::Size const cur_batch_size( multistage_vector_of_work[i].size() );
		if ( i == 1 || max_batch_size > cur_batch_size ) {
			max_batch_size = cur_batch_size;
		}
	}

	if ( max_batch_size == 0 ) {
		TR.Warning << "All work vectors in the multi-stage outer vector that was passed to the RosettaThreadManager have size zero!  Duly returning without doing anything." << std::endl;
		return;
	}

	platform::Size const actual_threads_to_request( std::min( max_batch_size, requested_thread_count ) );
	utility::vector1< utility::pointer::shared_ptr< utility::vector1< utility::thread::ReadWriteMutex > > > mutexes;
	utility::vector1< utility::vector1< bool > > jobs_completed( multistage_vector_of_work.size() );
	for ( platform::Size i(1), imax(multistage_vector_of_work.size()); i<=imax; ++i ) {
		debug_assert( i<=jobs_completed.size());
		platform::Size const curvect_size( multistage_vector_of_work[i].size() );
		utility::pointer::shared_ptr< utility::vector1< utility::thread::ReadWriteMutex > > mutvect( utility::pointer::make_shared< utility::vector1< utility::thread::ReadWriteMutex > >( curvect_size ) );
		mutexes.push_back( mutvect );
		jobs_completed[i].resize( curvect_size, false );
	}
	debug_assert(mutexes.size() == multistage_vector_of_work.size());
	debug_assert(jobs_completed.size() == multistage_vector_of_work.size());

	//Used to allow threads to block until all threads finish a grouped set of tasks:
	std::mutex barrier_mutex;
	utility::vector1< platform::Size > barrier_threadcount( multistage_vector_of_work.size() - 1, 0 );
	std::condition_variable barrier_cv;

	RosettaThreadAssignmentInfoOP thread_assignment2( thread_assignment == nullptr ? utility::pointer::make_shared< RosettaThreadAssignmentInfo >(RosettaThreadRequestOriginatingLevel::UNKNOWN) : thread_assignment );

	RosettaThreadFunctionOP fxn( utility::pointer::make_shared< RosettaThreadFunction > ( std::bind( &RosettaThreadManager::multistage_work_vector_thread_function, this, std::cref( multistage_vector_of_work ), std::ref( mutexes ), std::ref( jobs_completed ),
		std::ref( barrier_mutex ), std::ref( barrier_threadcount ), std::ref( barrier_cv ), thread_assignment2 ) ) );

	run_function_in_threads( fxn, actual_threads_to_request, RosettaThreadManagerAdvancedAPIKey(), thread_assignment2 );
#else
	//In non-threaded builds, just do the work.
	runtime_assert_string_msg( requested_thread_count == 0 || requested_thread_count == 1, "Error in RosettaThreadManager::do_multistage_work_vector_in_threads(): In non-threaded builds, only one thread may be requested.  Compile Rosetta with the \"extras=cxx11thread\" option to enable multi-threading." );
	for ( platform::Size i(1), imax(multistage_vector_of_work.size()); i<=imax; ++i ) {
		for ( platform::Size j(1), jmax(multistage_vector_of_work[i].size()); j<=jmax; ++j ) {
			(*(multistage_vector_of_work[i][j]))(); //Do jth piece of work in the ith round of work.
		}
	}
#endif
}

/// @brief ADVANCED API THAT SHOULD NOT BE USED IN MOST CIRCUMSTANCES.  Given a function that was bundled with its
/// arguments with std::bind, run it in many threads.  This calls RosettaThreadPool::run_function_in_threads for
/// the already-running thread pool.  If the thread pool has not been created, it first creates it by calling
/// create_thread_pool().  IF YOU DECIDE TO USE THE ADVANCED API, YOU MUST:
/// 1. Pass this function a RosettaThreadManagerAdvancedAPIKey from the calling context.  Since
/// the RosettaThreadManagerAdvancedAPIKey class has a private constructor, it can only be created in
/// whitelisted contexts in its friend list, which means that you must:
/// 2. Add the class that calls this advanced API to the friend list for the RosettaThreadManagerAdvancedAPIKey
/// class.  Since this will trigger breakage of the central_class_modification regression test, you must finally:
/// 3. Justify to the developer community why you must call this interface and not the safer, basic interface (do_work_vector_in_threads)
/// in both the comments in RosettaThreadManagerAdvancedAPIKey's friend list, the comments in the calling class, AND in your
/// pull request description.  Andrew Leaver-Fay and Vikram K. Mulligan will both scrutinize this closely.  It is highly recommended
/// that before using the run_function_in_threads() function, you first contact Andrew or Vikram and discuss whether it is possible
/// to do what you want to do using the basic API (the do_work_vector_in_threads() function).
///
/// @details The function is assigned to as many threads as the RosettaThreadPool decides to assign it to, always
/// including the thread from which the request originates.  It is guaranteed to run in 1 <= actual_thread_count <=
/// requested_thread_count threads.  After assigning the function to up to (requsted_thread_count - 1) other threads,
/// the function executes in the current thread, then the current thread blocks until the assigned threads report that
/// they are idle.  All of this is handled by the RosettaThreadPool class (or its derived classes, which may have)
/// different logic for assigning thread requests to threads).
///
/// @note Optionally, a RosettaThreadAssignmentInfo object can be passed in.  If provided, it will be populated with
/// the number of threads requested, the number actually assigned, the indices of the assigned threads, and a map of
/// system thread ID to Rosetta thread index.  The same owning pointer may optionally be provided to the function to
/// execute by the calling function if the function to execute requires access to this information.  Note also that the
/// function passed in is responsible for ensuring that it is able to carry out a large block of work, alone or concurrently
/// with many copies of itself in parallel threads, in a threadsafe manner.  Finally, note that this function requires a
/// RosettaThreadManagerAdvancedAPIKey, which can only be instantiated by friend classes in the whitelist in the
/// RosettaThreadManagerAdvancedAPIKey class definition.  This ensures that only select classes can access the advanced
/// RosettaThreadManager API.
void
RosettaThreadManager::run_function_in_threads(
	RosettaThreadFunctionOP function_to_execute,
	platform::Size const requested_thread_count,
	RosettaThreadManagerAdvancedAPIKey const &,
#ifdef MULTI_THREADED
	RosettaThreadAssignmentInfoOP thread_assignment /*= nullptr*/
#else
	RosettaThreadAssignmentInfoOP /*= nullptr*/
#endif
) {
#ifdef MULTI_THREADED
	{
		std::lock_guard< std::mutex > lock( thread_pool_mutex_ );
		if ( thread_pool_ == nullptr ) {
			create_thread_pool();
		}
	} //End of lock scope.
	thread_pool_->run_function_in_threads( function_to_execute, requested_thread_count, thread_assignment );
#else
	// In non-threaded builds, just run the function.
	runtime_assert_string_msg( requested_thread_count == 0 || requested_thread_count == 1, "Error in RosettaThreadManager::run_function_in_threads(): In non-threaded builds, only one thread may be requested.  Compile Rosetta with the \"extras=cxx11thread\" option to enable multi-threading." );
	(*function_to_execute)();
#endif
}

/// @brief Get the Rosetta thread index.
platform::Size
RosettaThreadManager::get_rosetta_thread_index() const {
#ifdef MULTI_THREADED
	{
		std::lock_guard< std::mutex > lock( thread_pool_mutex_ );
		if ( thread_pool_ == nullptr ) {
			return 0;
		}
	} //End of lock scope.
	if ( thread_id_to_rosetta_thread_index_.count( std::this_thread::get_id() ) ) {
		return thread_id_to_rosetta_thread_index_.at( std::this_thread::get_id() );
	} else {
		//Note: the following can't be a Rosetta tracer, since the Rosetta tracer asks for the Rosetta thread index:
		std::cerr << "Rosetta thread index was requested, but this is neither a thread that Rosetta launched, nor the thread from which Rosetta launched threads!  This can happen if external code (e.g. PyRosetta) has already launched threads." << std::endl;
		return 0;
	}
#else
	//In non-threaded builds, this always is thread 0.
	return 0;
#endif
}


// Private fxns /////////////////////////////////////////////////////

#ifdef MULTI_THREADED

/// @brief The function that is passed by do_work_vector_in_threads() to run_function_in_threads() to run in parallel,
/// to execute a vector of work in a threadsafe manner.
void
RosettaThreadManager::work_vector_thread_function(
	utility::vector1< RosettaThreadFunctionOP > const & vector_of_work,
	utility::vector1< utility::thread::ReadWriteMutex > & job_mutexes,
	utility::vector1< bool > & jobs_completed
) const {
	platform::Size jobcount(0);
	for ( platform::Size i(1), imax(vector_of_work.size()); i<=imax; ++i ) {
		{ //First, check with a read lock to see if this job is done.
			utility::thread::ReadLockGuard readlock( job_mutexes[i] );
			if ( jobs_completed[i] ) continue;
		}
		{ //Check again with a write lock, and set the job to completed if it isn't already completed.
			utility::thread::WriteLockGuard writelock( job_mutexes[i] );
			if ( jobs_completed[i] ) continue;
			jobs_completed[i] = true;
		}
		//If we reach here, we have a free hand to carry out the job.  (We don't even need a read lock, since the bool is flipped).
		(*vector_of_work[i])(); //Do the ith piece of work.
		if ( TR.visible() ) ++jobcount;
	}
	TR << "Thread " << get_rosetta_thread_index() << " completed " << jobcount << " of " << vector_of_work.size() << " work units." << std::endl; //Switch to debug output later.
}

/// @brief The function that is passed by do_multistage_work_vector_in_threads() to run_function_in_threads() to run in parallel,
/// to execute a vector of work in a threadsafe manner.
void
RosettaThreadManager::multistage_work_vector_thread_function(
	utility::vector1< utility::vector1< RosettaThreadFunctionOP > > const & multistage_vector_of_work,
	utility::vector1< utility::pointer::shared_ptr< utility::vector1< utility::thread::ReadWriteMutex > > > & multistage_job_mutexes,
	utility::vector1< utility::vector1< bool > > & multistage_jobs_completed,
	std::mutex & barrier_mutex,
	utility::vector1< platform::Size > & barrier_threadcount,
	std::condition_variable & barrier_cv,
	RosettaThreadAssignmentInfoOP thread_assignment
) const {
	platform::Size jobcount(0);

	platform::Size nthreads( thread_assignment->get_assigned_total_thread_count() );

	//Loop through the blocks of work that must be done sequentially:
	for ( platform::Size iblock(1), iblockmax( multistage_vector_of_work.size() ); iblock <= iblockmax; ++iblock ) {
		utility::vector1< RosettaThreadFunctionOP > const & vector_of_work( multistage_vector_of_work[iblock] );
		utility::vector1< utility::thread::ReadWriteMutex > & job_mutexes( *(multistage_job_mutexes[iblock]) );
		utility::vector1< bool > & jobs_completed( multistage_jobs_completed[iblock] );

		debug_assert( vector_of_work.size() == job_mutexes.size() );
		debug_assert( vector_of_work.size() == jobs_completed.size() );

		for ( platform::Size i(1), imax(vector_of_work.size()); i<=imax; ++i ) {
			{ //First, check with a read lock to see if this job is done.
				utility::thread::ReadLockGuard readlock( job_mutexes[i] );
				if ( jobs_completed[i] ) continue;
			}
			{ //Check again with a write lock, and set the job to completed if it isn't already completed.
				utility::thread::WriteLockGuard writelock( job_mutexes[i] );
				if ( jobs_completed[i] ) continue;
				jobs_completed[i] = true;
			}
			//If we reach here, we have a free hand to carry out the job.  (We don't even need a read lock, since the bool is flipped).
			(*vector_of_work[i])(); //Do the ith piece of work.
			if ( TR.visible() ) ++jobcount;
		}
		TR << "In work block " << iblock << " of " << iblockmax << ", thread " << get_rosetta_thread_index() << " completed " << jobcount << " of " << vector_of_work.size() << " work units." << std::endl; //Switch to debug output later.
		if ( iblock < iblockmax ) {
			bool triggering_thread(false);
			{ //Scope for mutex
				std::lock_guard< std::mutex > lock( barrier_mutex );
				++barrier_threadcount[iblock];
				if ( barrier_threadcount[iblock] == nthreads ) {
					//TR << "Thread " << get_rosetta_thread_index() << " is the last one to the barrier after round " << iblock << "." << std::endl; //DELETE ME
					barrier_cv.notify_all();
					triggering_thread = true;
				}
			} //End scope for mutex.
			if ( !triggering_thread ) {
				//TR << "Thread " << get_rosetta_thread_index() << " is waiting at the barrier after round " << iblock << "." << std::endl; //DELETE ME.
				std::unique_lock< std::mutex > uniquelock( barrier_mutex );
				barrier_cv.wait( uniquelock, [&barrier_threadcount, iblock, nthreads]{ return ( barrier_threadcount[iblock] == nthreads ); } );
			}
			//TR << "*Click*.  Thread " << get_rosetta_thread_index() << " is past the barrier after round " << iblock << "." << std::endl; //DELETE ME.
		}
	}
}

#endif

} //thread_manager
} //basic
