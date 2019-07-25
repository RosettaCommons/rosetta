// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadPool.cc
/// @brief A container for a pool of threads.  Threads idle continuously until loaded with a
/// function to execute.  A subset of threads then execute the function synchronously and then
/// return to the idle state.
/// @note The basic verison of this class assigns functions to threads on a first-come, first-served
/// basis.  If I request 8 threads and 8 threads are idle, I get all 8 threads.  If the next request
/// also asks for 8 threads while I'm still using mine, it only gets the thread that made the request.
/// Derived classes can implement different behaviour for assigning work to threads.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef MULTI_THREADED

#include <basic/thread_manager/RosettaThreadPool.hh>
#include <basic/thread_manager/RosettaThreadAssignmentInfo.hh>
#include <basic/thread_manager/RosettaThread.hh>
#include <basic/Tracer.hh>

#include <condition_variable>

static basic::Tracer TR( "basic.thread_manager.RosettaThreadPool" );


namespace basic {
namespace thread_manager {

/// @brief Constructor requiring the number of threads to launch.
/// @note The pool holds total_threads-1 threads.  A RosettaThreadPoolInstantiationKey must be passed to prove that
/// this is being instantiated in a context that is allowed to instantiate the RosettaThreadPool (i.e. the context of
/// the RosettaThreadManager).
RosettaThreadPool::RosettaThreadPool(
	platform::Size const total_threads,
	RosettaThreadPoolInstantiationKey const &
) :
	utility::pointer::ReferenceCount()
{
	runtime_assert_string_msg( total_threads > 0, "Error in RosettaThreadPool constructor: The number of threads must be greater than or equal to 1." );
	threads_.reserve( total_threads - 1 );
	for ( platform::Size i(1); i<total_threads; ++i ) {
		threads_.push_back( utility::pointer::make_shared< RosettaThread >( i, RosettaThreadInstantiationKey() ) );
	}
	TR << "Launched " << total_threads - 1 << " new threads." << std::endl;
	debug_assert( threads_.size() == total_threads - 1 );
}

/// @brief Destructor.
RosettaThreadPool::~RosettaThreadPool(){
	force_terminate_threads();
}

/// @brief Given a function that was bundled with its arguments with std::bind, run it in many threads.
/// @details The function is assigned to as many threads as are available, including the thread from which the request
/// originates.  It is guaranteed to run in 1 <= actual_thread_count <= requested_thread_count threads.  After assigning
/// the function to up to (requsted_thread_count - 1) other threads, the function executes in the current thread, then
/// the current thread blocks until the assigned threads report that they are idle.
/// @note Optionally, a RosettaThreadAssignmentInfo object can be passed in.  If provided, it will be populated with
/// the number of threads requested, the number actually assigned, the indices of the assigned threads, etc.  The same
/// owning pointer may optionally be provided to the function to execute by the calling function if the function to
/// execute requires access to this information.
void
RosettaThreadPool::run_function_in_threads(
	RosettaThreadFunctionOP function_to_execute,
	Size const requested_thread_count,
	RosettaThreadAssignmentInfoOP thread_assignment /*=nullptr*/
) {
	runtime_assert_string_msg( function_to_execute != nullptr, "Error in RosettaThreadPool::run_function_in_threads():  A null pointer was passed to this function!" );
	runtime_assert_string_msg( requested_thread_count > 0, "Error in RosettaThreadPool::run_function_in_threads():  Zero threads were requested!");
	if ( requested_thread_count > threads_.size() + 1 ) {
		TR.Warning << "Warning: " << requested_thread_count << " threads were requested, but only " << threads_.size() + 1 << " exist in the RosettaThreadPool." << std::endl;
	}

	utility::vector1< platform::Size > assigned_threads;
	platform::Size assigned_thread_count( 0 );

	//Variables to use to determine when all jobs are done:
	std::mutex job_completion_mutex;
	std::unique_lock< std::mutex > job_completion_lock( job_completion_mutex );
	platform::Size jobs_finished(0);
	std::condition_variable cv;

	//Assign work to other threads.  Note that if the request is coming from one of the threads in the pool,
	//that thread will be considered "not idle", so work won't be assigned to it in this step.  Below, the
	//function will *also* run in this thread.
	if ( requested_thread_count > 1 ) {
		//We lock the thread assignment mutex while assigning to threads, to prevent other threads from trying to assign to the same threads.
		std::lock_guard< std::mutex > lock( working_thread_list_mutex_ );
		for ( platform::Size i(1), imax(threads_.size()); i<=imax; ++i ) {
			if ( threads_[i]->is_available_for_new_work() && threads_[i]->is_idle() ) {
				threads_[i]->set_forced_idle(true);
				runtime_assert( threads_[i]->set_function( function_to_execute, job_completion_mutex, cv, jobs_finished ) ); //Should always succeed.
				++assigned_thread_count;
				assigned_threads.push_back( i );
				if ( assigned_thread_count == requested_thread_count - 1 ) break;
			}
		}

		//Store the information about assignments
		if ( thread_assignment != nullptr ) {
			debug_assert( assigned_thread_count == assigned_threads.size() );
			thread_assignment->set_assigned_child_threads( assigned_threads );
			thread_assignment->set_requested_thread_count( requested_thread_count );
		}

		//Release all of the threads and let them run
		for ( platform::Size const & thread_index : assigned_threads ) {
			threads_[thread_index]->set_forced_idle(false);
		}
	} //End of mutex lock scope

	//Also run the function in this thread (the requesting thread):
	(*function_to_execute)();

	//Now wait until all the other threads have finished.
	if ( assigned_thread_count > 0 ) {

		//Wait for all threads to finish
		if ( jobs_finished != assigned_thread_count ) {
			cv.wait( job_completion_lock, [&jobs_finished, assigned_thread_count]{ return jobs_finished == assigned_thread_count; } );
		}

		//Once all threads have finished, mark them all as available for new work.
		{
			std::lock_guard< std::mutex > lock( working_thread_list_mutex_ ); //Don't allow work assignments.
			for ( platform::Size const & thread_index : assigned_threads ) {
				threads_[thread_index]->set_available_for_new_work();
			}
		} //End of lock scope.
	}
} //run_function_in_threads

/// @brief Force all threads to terminate.  Called by the destructor of the RosettaThreadManager class, and
/// by the destructor of this class.
/// @details At the end of this operation, the threads_ object is empty.
void
RosettaThreadPool::force_terminate_threads() {
	std::lock_guard< std::mutex > lock( working_thread_list_mutex_ ); //Obtain a lock to delete threads.
	if ( threads_.size() > 0 ) {
		//TR << "Initiating termination of " << threads_.size() << " threads." << std::endl;
		threads_.clear();
	}
	//TR << "All threads terminated." << std::endl;
}

/// @brief Construct a map of (system thread id)->(Rosetta thread index), and return the object.
/// @details Intended to be called only by the RosettaThreadManager.  If you want the current thread ID, use
/// basic::thread_manager::RosettaThreadManager::get_instance()->get_rosetta_thread_index().
std::map< std::thread::id, platform::Size >
RosettaThreadPool::get_thread_id_to_rosetta_thread_index_map() const {
	std::map< std::thread::id, platform::Size > return_map;
	return_map[ std::this_thread::get_id() ] = 0;
	std::lock_guard< std::mutex > lock( working_thread_list_mutex_ ); //Obtain a lock to query threads, to be on the safe side.
	for ( platform::Size i(1), imax( threads_.size() ); i<=imax; ++i ) {
		return_map[ threads_[i]->get_system_thread_id() ] = threads_[i]->get_thread_index();
	}
	debug_assert( return_map.size() == threads_.size() + 1 );
	return return_map;
}


} //thread_manager
} //basic

#endif //MULTI_THREADED
