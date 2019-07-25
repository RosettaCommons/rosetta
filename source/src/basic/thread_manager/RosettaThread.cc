// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThread.cc
/// @brief A container for a thread in a RosettaThreadPool.  The thread idles continuously until
/// loaded with a function to execute.  It then executes the function and returns to the idle state.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef MULTI_THREADED

#include <basic/thread_manager/RosettaThread.hh>
#include <utility/sys_util.hh>
#include <basic/Tracer.hh>
#include <basic/random/init_random_generator.hh>

static basic::Tracer TR( "basic.thread_manager.RosettaThread" );

namespace basic {
namespace thread_manager {

/// @brief RosettaThreadPool constructor.  The thread index *must* be passed.
RosettaThread::RosettaThread(
	platform::Size const thread_index,
	RosettaThreadInstantiationKey const &
) :
	utility::pointer::ReferenceCount(),
	thread_index_( thread_index )
{
	hold_thread_ = false;
	is_available_for_new_work_ = true;
	terminate_now_ = false;
	runtime_assert_string_msg( thread_index_ != 0, "Error in RosettaThread constructor:  The thread index cannot be zero!" );
	this_thread_ = std::thread( &RosettaThread::thread_function, this ); //Launch a thread.
}

/// @brief Destructor.  Calls force_thread_termination().
RosettaThread::~RosettaThread(){
	force_thread_termination();
}

//////////////// PUBLIC FUNCTIONS ////////////////

/// @brief Set the function to run on this thread.
/// @details Requires a function that returns void (bundled with its arguments with std::bind).
bool
RosettaThread::set_function(
	RosettaThreadFunctionOP function_to_execute,
	std::mutex & job_completion_mutex,
	std::condition_variable & job_completion_condition_variable,
	platform::Size & job_completion_count
) {
	runtime_assert_string_msg( function_to_execute != nullptr, "Error in RosettaThread::set_function(): The function pointer cannot be null." );
	{
		std::lock_guard< std::mutex > locker( this_thread_function_mutex_ );
		if ( function_to_execute_ != nullptr ) return false; //Shouldn't be possible.  Might make this an error.
		is_available_for_new_work_ = false; //Atomic operation.
		function_to_execute_ = function_to_execute;
		job_completion_mutex_ = &job_completion_mutex;
		job_completion_condition_variable_ = &job_completion_condition_variable;
		job_completion_count_ = &job_completion_count;
	}
	cv_.notify_one();
	return true;
}

/// @brief Set the idle state of this thread.  If set to true and a function is already running, it continues
/// to run, but no new function will run until this is set to false.
void
RosettaThread::set_forced_idle(
	bool const setting
) {
	{
		std::lock_guard< std::mutex > lock( this_thread_function_mutex_ );
		hold_thread_ = setting; //Atomic operation.
	}
	if ( setting == false ) cv_.notify_one();
}

/// @brief Is this thread idle?
bool
RosettaThread::is_idle() const {
	std::lock_guard< std::mutex > lock( this_thread_function_mutex_ );
	return function_to_execute_ == nullptr;
}

//////////////// PRIVATE FUNCTIONS ////////////////

/// @brief The function that this thread performs, which idles until the function object is populated.
void
RosettaThread::thread_function() {
	TR << "Launching thread " << thread_index_ << "." << std::endl;
	basic::random::init_random_number_generators(thread_index_);

	do {
		std::unique_lock< std::mutex > uniquelock( this_thread_function_mutex_ );
		if ( !( terminate_now_.load() || (!(hold_thread_.load()) && function_to_execute_ != nullptr) ) ) {
			cv_.wait( uniquelock, [this]{ return ( terminate_now_.load() || (!(hold_thread_.load()) && function_to_execute_ != nullptr) ); } );
		}
		//After waiting, the mutex is locked.

		if ( terminate_now_.load() ) { //Signal to spin down.
			break;
		}
		if ( !(hold_thread_.load()) /*Thread not externally set to idle.*/ && function_to_execute_ != nullptr /*We have a function to execute.*/ ) {
			(*function_to_execute_)(); //Run the function to execute in this thread.
			function_to_execute_.reset(); //Delete the function (or at least clear this owning pointer, if there are other owners).
			debug_assert( job_completion_mutex_ != nullptr );
			debug_assert( job_completion_condition_variable_ != nullptr );
			debug_assert( job_completion_count_ != nullptr );
			{
				std::lock_guard< std::mutex > lock( *job_completion_mutex_ );
				++(*job_completion_count_);
				job_completion_condition_variable_->notify_one(); //Moving this inside the scope of the lock to try to avoid a rare race condition.
			}

			job_completion_mutex_ = nullptr;
			job_completion_count_ = nullptr;
			job_completion_condition_variable_ = nullptr;

			debug_assert( function_to_execute_ == nullptr ); //Should be true, now.
			debug_assert( job_completion_mutex_ == nullptr );
			debug_assert( job_completion_condition_variable_ == nullptr );
			debug_assert( job_completion_count_ == nullptr );
		}
	} while(true); //Loop forever
} //thread_function()

/// @brief Terminates the std::thread object held by this object.
/// @details Warning!  This should only be called by this object's destructor.  Renders the object unusable afterwards.
void
RosettaThread::force_thread_termination() {
	{
		std::lock_guard< std::mutex > lock(this_thread_function_mutex_);
		if ( !is_idle_already_locked() ) {
			utility_exit_with_message( "RosettaThread::force_thread_termination():  Thread " + std::to_string(thread_index_) + " is not idle!  Cannot halt thread." );
		} else {
			//TR << "Terminating thread " << thread_index_ << " (which is idle)." << std::endl;
			terminate_now_ = true;
		}
	}
	cv_.notify_one();
	this_thread_.join();
}

} //thread_manager
} //basic

#endif //MULTI_THREADED
