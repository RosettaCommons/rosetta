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
#ifdef MULTI_THREADED
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/CitationManager.hh>
#endif

// C++ headers
#include <string>
#include <cmath>

// Construct tracer.
static basic::Tracer TR( "basic.thread_manager.RosettaThreadManager" );

namespace basic {
namespace thread_manager {

#define MAX_THREAD_INDEX_WARNINGS 8 //The maximum number of warnings about threads not managed by the RosettaThreadManager that the user will see.

// Private methods ////////////////////////////////////////////////////////////

/// @brief Empty constructor.
RosettaThreadManager::RosettaThreadManager()
#ifdef MULTI_THREADED
:
	nthreads_( basic::options::option[ basic::options::OptionKeys::multithreading::total_threads ]() ),
	thread_pool_was_launched_(false),
	warning_counter_(0)
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
/// @details Accesses the global options system to determine the number of threads to launch.  Also, registers the
/// RosettaThreadManager with the CitationManager if launching threads.
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
		register_thread_manager_with_citation_manager();
	} else {
		TR.Warning << "WARNING WARNING WARNING!  RosettaThreadManager::create_thread_pool() was called, but the thread pool has already been created!" << std::endl;
	}

	RosettaThreadManagerInitializationTracker::get_instance()->mark_thread_manager_as_initialized();
	debug_assert( thread_pool_was_launched_ == false );
	thread_pool_was_launched_ = true;
}

/// @brief Adds citation information for the RosettaThreadManager to the CitationManager.
/// @details The thread_pool_mutex_ must be locked before calling this function!
void
RosettaThreadManager::register_thread_manager_with_citation_manager() const {
	using namespace basic::citation_manager;
	CitationManager::get_instance()->add_unpublished_modules(
		utility::vector1< UnpublishedModuleInfoCOP > {
		utility::pointer::make_shared< UnpublishedModuleInfo >(
		"RosettaThreadManager", CitedModuleType::Singleton,
		"Vikram K. Mulligan",
		"Systems Biology, Center for Computational Biology, Flatiron Institute",
		"vmulligan@flatironinstitute.org"
		)
		}
	);
}

/// @brief Trigger launch of threads.
/// @details Does nothing if threads already launched.
void
RosettaThreadManager::launch_threads() {
	if ( thread_pool_was_launched_ == false ) {
		std::lock_guard< std::mutex > lock( thread_pool_mutex_ );
		if ( thread_pool_ == nullptr ) {
			create_thread_pool();
		}
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
	utility::vector1< RosettaThreadFunction > const & vector_of_work,
	platform::Size const requested_thread_count,
#ifdef MULTI_THREADED
	RosettaThreadAssignmentInfo & thread_assignment
#else
	RosettaThreadAssignmentInfo &
#endif
) {
	if ( vector_of_work.size() == 0 ) {
		TR.Warning << "A work vector of size zero was passed to the RosettaThreadManager!  Duly returning without doing anything." << std::endl;
		return;
	}
#ifdef MULTI_THREADED
	platform::Size const nonzero_requested_thread_count( requested_thread_count == 0 ? nthreads_ : requested_thread_count );
	platform::Size const actual_threads_to_request( std::min( vector_of_work.size(), nonzero_requested_thread_count ) );
	utility::vector1< std::mutex > mutexes( vector_of_work.size() );
	utility::vector1< AtomicBoolContainer > jobs_completed( vector_of_work.size() ); // Initialized to false automatically.

	RosettaThreadFunction fxn(
		std::bind(
		&RosettaThreadManager::work_vector_thread_function, this,
		std::cref( vector_of_work ), std::ref( mutexes ), std::ref( jobs_completed ), std::cref( thread_assignment )
		)
	);
	run_function_in_threads( fxn, actual_threads_to_request, RosettaThreadManagerAdvancedAPIKey(), thread_assignment );
#else
	//In non-threaded builds, just iterate through and do all of the work.
	runtime_assert_string_msg( requested_thread_count == 0 || requested_thread_count == 1, "Error in RosettaThreadManager::do_work_vector_in_threads(): In non-threaded builds, only one thread may be requested.  Compile Rosetta with the \"extras=cxx11thread\" option to enable multi-threading." );
	for ( platform::Size i(1), imax(vector_of_work.size()); i<=imax; ++i ) {
		(vector_of_work[i])(); //Do the ith piece of work.
	}
#endif
}

/// @brief VARIANT BASIC API THAT SHOULD BE USED FOR WORK VECTORS OF NEAR-EQUAL SIZED CHUNKS WHERE THE CHUNKS ARE SMALL.  Given a vector of
/// functions that were bundled with their arguments with std::bind, each of which can be executed in any order and each of which is safe to execute
/// in parallel with any other, run all of these in threads.
/// @details The bundled functions should be atomistic pieces of work.  They should be bundled with their arguments with
/// std::bind, and the arguments should include the place to store output (i.e. they should return void).  These functions
/// should not handle any mutexes themselves, but should ensure that they are operating only on memory locations that no
/// other functions in the vector are operating on.
/// @note Under the hood, this sets up no mutexes, instead giving each thread a staggered subset of the work in the vector.  It calls
/// run_function_in_threads() to do the work.  The work is done concurrently in 1 <= actual count <= min( requested thread count, total
/// thread count ) threads.  The function blocks until all threads have finished their work, which means that the individual work units
/// should be small, that the longest-running work unit should be short compared to the total runtime, and that the number of work units
/// should be much greater than the number of threads requested.
/// This function works best for cases in which it is known that most of the work in the vector is of equal size (i.e. load-balancing is
/// unlikely to be an issue), and where the overhead of locking mutexes for each job is likely to be comparable in size to the cost of a
/// job (so we want to avoid this overhead).
void
RosettaThreadManager::do_work_vector_in_threads_no_locking(
	utility::vector1< RosettaThreadFunction > const & vector_of_work,
	platform::Size const requested_thread_count,
#ifdef MULTI_THREADED
	RosettaThreadAssignmentInfo & thread_assignment
#else
	RosettaThreadAssignmentInfo &
#endif
) {
#ifdef MULTI_THREADED
	platform::Size const nonzero_requested_thread_count( requested_thread_count == 0 ? nthreads_ : requested_thread_count );
	platform::Size const actual_threads_to_request( std::min( vector_of_work.size(), nonzero_requested_thread_count ) );
	RosettaThreadFunction fxn(
		std::bind(
		&RosettaThreadManager::work_vector_thread_function_no_locking, this,
		std::cref( vector_of_work ), std::cref( thread_assignment )
		)
	);
	run_function_in_threads( fxn, actual_threads_to_request, RosettaThreadManagerAdvancedAPIKey(), thread_assignment );
#else
	//In non-threaded builds, just iterate through and do all of the work.
	runtime_assert_string_msg( requested_thread_count == 0 || requested_thread_count == 1, "Error in RosettaThreadManager::do_work_vector_in_threads_no_locking(): In non-threaded builds, only one thread may be requested.  Compile Rosetta with the \"extras=cxx11thread\" option to enable multi-threading." );
	for ( platform::Size i(1), imax(vector_of_work.size()); i<=imax; ++i ) {
		(vector_of_work[i])(); //Do the ith piece of work.
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
	utility::vector1< utility::vector1< RosettaThreadFunction > > const & multistage_vector_of_work,
	platform::Size const requested_thread_count,
#ifdef MULTI_THREADED
	RosettaThreadAssignmentInfo & thread_assignment
#else
	RosettaThreadAssignmentInfo &
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
	std::deque< std::deque< std::mutex > > mutexes;
	std::deque< std::deque< AtomicBoolContainer > > jobs_completed;
	for ( platform::Size i(0), imax(multistage_vector_of_work.size()); i<imax; ++i ) {
		platform::Size const curvect_size( multistage_vector_of_work[i+1].size() );
		mutexes.emplace_back( curvect_size );
		jobs_completed.emplace_back( curvect_size ); //Initialized to false automatically.
	}
	debug_assert(mutexes.size() == multistage_vector_of_work.size());
	debug_assert(jobs_completed.size() == multistage_vector_of_work.size());

	//Used to allow threads to block until all threads finish a grouped set of tasks:
	std::mutex barrier_mutex;
	utility::vector1< platform::Size > barrier_threadcount( multistage_vector_of_work.size() - 1, 0 );
	std::condition_variable barrier_cv;

	RosettaThreadFunction fxn(
		std::bind(
		&RosettaThreadManager::multistage_work_vector_thread_function, this,
		std::cref( multistage_vector_of_work ), std::ref( mutexes ), std::ref( jobs_completed ),
		std::ref( barrier_mutex ), std::ref( barrier_threadcount ), std::ref( barrier_cv ), std::ref(thread_assignment)
		)
	);

	run_function_in_threads( fxn, actual_threads_to_request, RosettaThreadManagerAdvancedAPIKey(), thread_assignment );
#else
	//In non-threaded builds, just do the work.
	runtime_assert_string_msg( requested_thread_count == 0 || requested_thread_count == 1, "Error in RosettaThreadManager::do_multistage_work_vector_in_threads(): In non-threaded builds, only one thread may be requested.  Compile Rosetta with the \"extras=cxx11thread\" option to enable multi-threading." );
	for ( platform::Size i(1), imax(multistage_vector_of_work.size()); i<=imax; ++i ) {
		for ( platform::Size j(1), jmax(multistage_vector_of_work[i].size()); j<=jmax; ++j ) {
			((multistage_vector_of_work[i][j]))(); //Do jth piece of work in the ith round of work.
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
/// @note A RosettaThreadAssignmentInfo object should be passed in.  It will be populated with
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
	RosettaThreadFunction & function_to_execute,
	platform::Size const requested_thread_count,
	RosettaThreadManagerAdvancedAPIKey const &,
#ifdef MULTI_THREADED
	RosettaThreadAssignmentInfo & thread_assignment
#else
	RosettaThreadAssignmentInfo &
#endif
) {
#ifdef MULTI_THREADED
	if ( thread_pool_was_launched_ == false ) {
		std::lock_guard< std::mutex > lock( thread_pool_mutex_ );
		if ( thread_pool_ == nullptr ) {
			create_thread_pool();
		}
	} //End of lock scope.
	thread_pool_->run_function_in_threads( &function_to_execute, requested_thread_count, thread_assignment );
#else
	// In non-threaded builds, just run the function.
	runtime_assert_string_msg( requested_thread_count == 0 || requested_thread_count == 1, "Error in RosettaThreadManager::run_function_in_threads(): In non-threaded builds, only one thread may be requested.  Compile Rosetta with the \"extras=cxx11thread\" option to enable multi-threading." );
	(function_to_execute)();
#endif
}

/// @brief Get the Rosetta thread index.
platform::Size
RosettaThreadManager::get_rosetta_thread_index() const {
#ifdef MULTI_THREADED
	if ( thread_pool_was_launched_ == false ) {
		return 0;
	} //End of lock scope.
	if ( thread_id_to_rosetta_thread_index_.count( std::this_thread::get_id() ) ) {
		return thread_id_to_rosetta_thread_index_.at( std::this_thread::get_id() );
	} else {
		//Note: the following can't be a Rosetta tracer, since the Rosetta tracer asks for the Rosetta thread index:
		if ( warning_counter_ < MAX_THREAD_INDEX_WARNINGS ) {
			std::lock_guard<std::mutex> warning_lock( warning_counter_mutex_ ); //Lock the mutex for the scope of the if statement
			if ( warning_counter_ < MAX_THREAD_INDEX_WARNINGS ) {
				++warning_counter_;
				std::cerr << "Rosetta thread index was requested, but this is neither a thread that Rosetta launched, nor the thread from which Rosetta launched threads!  This can happen if external code (e.g. PyRosetta) has already launched threads.";
				if ( warning_counter_ >= MAX_THREAD_INDEX_WARNINGS ) {
					std::cerr << "  Silencing further instances of this warning.";
				}
				std::cerr << std::endl;
			}
		}
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
	utility::vector1< RosettaThreadFunction > const & vector_of_work,
	utility::vector1< std::mutex > & job_mutexes,
	utility::vector1< AtomicBoolContainer > & jobs_completed,
	RosettaThreadAssignmentInfo const & thread_assignments
) const {
	platform::Size jobcount(0);
	bool const debug_vis( TR.Debug.visible() );

	// The following is a little efficiency trick: each thread looks through the whole work vector for work to do, but we start
	// threads at staggered positions throughout the work vector so that they're not likely to be encountering blocks of work that
	// other threads have done.
	platform::Size const nthreads( thread_assignments.get_assigned_total_thread_count() );
	platform::Size const this_thread( thread_assignments.get_this_thread_index_in_assigned_set() );
	platform::Size const startpos( std::round( (this_thread - 1) * float( vector_of_work.size() ) / float( nthreads ) ) ); //Deliberately casting to float and not double.
	platform::Size const endpos( startpos == 0 ? vector_of_work.size() : startpos );

	platform::Size i(startpos);
	do {
		++i;
		if ( i > vector_of_work.size() ) i = 1;

		{ //First, check if this job is done.
			if ( jobs_completed[i].contained_bool_ == true ) continue;
		}
		{ //Check again with a lock, and set the job to completed if it isn't already completed.
			std::lock_guard< std::mutex > lock( job_mutexes[i] );
			if ( jobs_completed[i].contained_bool_ == true ) continue;
			jobs_completed[i].contained_bool_ = true;
		}
		//If we reach here, we have a free hand to carry out the job.  (We don't even need a read lock, since the bool is flipped).
		(vector_of_work[i])(); //Do the ith piece of work.
		if ( debug_vis ) ++jobcount;
	} while (i != endpos);
	if ( debug_vis ) {
		TR.Debug << "Thread " << get_rosetta_thread_index() << " completed " << jobcount << " of " << vector_of_work.size() << " work units." << std::endl;
	}
}

/// @brief The function that is passed by do_work_vector_in_threads_no_locking() to run_function_in_threads() to run in parallel,
/// to execute a vector of work in a threadsafe manner, without locking each task.
/// @details This version assigns every Nth piece of work to a given thread.  The assumption is that this will result in even load-balancing
/// without the overhead of locking.  This is true if the pieces of work are of roughly equal size.
void
RosettaThreadManager::work_vector_thread_function_no_locking(
	utility::vector1< RosettaThreadFunction > const & vector_of_work,
	RosettaThreadAssignmentInfo const & thread_assignments
) const {
	platform::Size jobcount(0);
	bool const debug_vis( TR.Debug.visible() );

	// The following is used to stagger the jobs done by each thread:
	platform::Size const nthreads( thread_assignments.get_assigned_total_thread_count() );
	platform::Size const this_thread( thread_assignments.get_this_thread_index_in_assigned_set() );

	for ( platform::Size i(this_thread), imax( vector_of_work.size() ); i<=imax; i+=nthreads ) { //Do every Nth piece of work.
		(vector_of_work[i])(); //Do the ith piece of work.
		if ( debug_vis ) ++jobcount;
	}
	if ( debug_vis ) {
		TR.Debug << "Thread " << get_rosetta_thread_index() << " completed " << jobcount << " of " << vector_of_work.size() << " work units." << std::endl;
	}
}

/// @brief The function that is passed by do_multistage_work_vector_in_threads() to run_function_in_threads() to run in parallel,
/// to execute a vector of work in a threadsafe manner.
void
RosettaThreadManager::multistage_work_vector_thread_function(
	utility::vector1< utility::vector1< RosettaThreadFunction > > const & multistage_vector_of_work,
	std::deque< std::deque< std::mutex > > & multistage_job_mutexes,
	std::deque< std::deque< AtomicBoolContainer > > & multistage_jobs_completed,
	std::mutex & barrier_mutex,
	utility::vector1< platform::Size > & barrier_threadcount,
	std::condition_variable & barrier_cv,
	RosettaThreadAssignmentInfo & thread_assignment
) const {
	platform::Size jobcount(0);

	// The following is a little efficiency trick: each thread looks through the whole work vector for work to do, but we start
	// threads at staggered positions throughout the work vector so that they're not likely to be encountering blocks of work that
	// other threads have done.);
	platform::Size const nthreads( thread_assignment.get_assigned_total_thread_count() );
	platform::Size const this_thread( thread_assignment.get_this_thread_index_in_assigned_set() );

	bool const debug_vis( TR.Debug.visible() );

	//Loop through the blocks of work that must be done sequentially:
	for ( platform::Size iblock(1), iblockmax( multistage_vector_of_work.size() ); iblock <= iblockmax; ++iblock ) {
		utility::vector1< RosettaThreadFunction > const & vector_of_work( multistage_vector_of_work[iblock] );
		std::deque< std::mutex > & job_mutexes( multistage_job_mutexes[iblock - 1] );
		std::deque< AtomicBoolContainer > & jobs_completed( multistage_jobs_completed[iblock - 1] );

		debug_assert( vector_of_work.size() == job_mutexes.size() );
		debug_assert( vector_of_work.size() == jobs_completed.size() );

		platform::Size const startpos( std::round( (this_thread - 1) * float( vector_of_work.size() ) / float( nthreads ) ) ); //Deliberately casting to float and not double.
		platform::Size const endpos( startpos == 0 ? vector_of_work.size() : startpos );

		platform::Size i(startpos);
		do {
			++i;
			if ( i > vector_of_work.size() ) {
				i = 1;
			}

			{ //First, check to see if this job is done.
				if ( jobs_completed[i-1].contained_bool_ == true ) continue;
			}
			{ //Check again with a lock, and set the job to completed if it isn't already completed.
				std::lock_guard< std::mutex > lock( job_mutexes[i-1] );
				if ( jobs_completed[i-1].contained_bool_ == true ) continue;
				jobs_completed[i-1].contained_bool_ = true;
			}
			//If we reach here, we have a free hand to carry out the job.  (We don't even need a read lock, since the bool is flipped).
			(vector_of_work[i])(); //Do the ith piece of work.
			if ( debug_vis ) ++jobcount;
		} while ( i != endpos );
		if ( debug_vis ) {
			TR.Debug << "In work block " << iblock << " of " << iblockmax << ", thread " << get_rosetta_thread_index() << " completed " << jobcount << " of " << vector_of_work.size() << " work units." << std::endl; //Switch to debug output later.
		}
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
