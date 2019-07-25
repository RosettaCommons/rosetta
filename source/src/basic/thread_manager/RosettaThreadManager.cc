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
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef MULTI_THREADED

#include <basic/thread_manager/RosettaThreadManager.hh>
#include <basic/thread_manager/RosettaThreadPool.hh>

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
{}

/// @brief Destructor.  Non-empty, since threads must be spun down.
RosettaThreadManager::~RosettaThreadManager() {
	if ( thread_pool_ != nullptr ) {
		//TR << "Forcing termination of threads that are still running." << std::endl;
		thread_pool_->force_terminate_threads();
	}
}

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
		platform::Size const nthreads( option[multithreading::total_threads] );
		TR << "Creating a thread pool of " << nthreads << " threads." << std::endl;
		thread_pool_ = utility::pointer::make_shared< RosettaThreadPool >( nthreads, RosettaThreadPoolInstantiationKey() );
		thread_id_to_rosetta_thread_index_ = thread_pool_->get_thread_id_to_rosetta_thread_index_map();
	} else {
		TR.Warning << "WARNING WARNING WARNING!  RosettaThreadManager::create_thread_pool() was called, but the thread pool has already been created!" << std::endl;
	}

	RosettaThreadManagerInitializationTracker::get_instance()->mark_thread_manager_as_initialized();
}

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
	RosettaThreadAssignmentInfoOP thread_assignment /*= nullptr*/
) {
	platform::Size const actual_threads_to_request( std::min( vector_of_work.size(), requested_thread_count ) );
	utility::vector1< utility::thread::ReadWriteMutex > mutexes( vector_of_work.size() );
	utility::vector1< bool > jobs_completed( vector_of_work.size(), false );

	RosettaThreadFunctionOP fxn( utility::pointer::make_shared< RosettaThreadFunction > ( std::bind( &RosettaThreadManager::work_vector_thread_function, this, std::cref( vector_of_work ), std::ref( mutexes ), std::ref( jobs_completed ) ) ) );
	run_function_in_threads( fxn, actual_threads_to_request, RosettaThreadManagerAdvancedAPIKey(), thread_assignment );
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
	RosettaThreadAssignmentInfoOP thread_assignment /*=nullptr*/
) {
	{
		std::lock_guard< std::mutex > lock( thread_pool_mutex_ );
		if ( thread_pool_ == nullptr ) {
			create_thread_pool();
		}
	} //End of lock scope.
	thread_pool_->run_function_in_threads( function_to_execute, requested_thread_count, thread_assignment );
}

/// @brief Get the Rosetta thread index.
platform::Size
RosettaThreadManager::get_rosetta_thread_index() const {
	return thread_id_to_rosetta_thread_index_.at( std::this_thread::get_id() );
}


// Private fxns /////////////////////////////////////////////////////

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

} //thread_manager
} //basic

#endif //MULTI_THREADED
