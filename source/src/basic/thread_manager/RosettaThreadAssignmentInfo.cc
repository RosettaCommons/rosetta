// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadAssignmentInfo.cc
/// @brief A class for storing information about the threads to which a function has been assigned by the RosettaThreadManager.
/// @details Accessor methods are fully threadsafe.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef MULTI_THREADED

#include <basic/thread_manager/RosettaThreadAssignmentInfo.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "basic.thread_manager.RosettaThreadAssignmentInfo" );


namespace basic {
namespace thread_manager {

/// @brief Default constructor.
RosettaThreadAssignmentInfo::RosettaThreadAssignmentInfo(
	RosettaThreadRequestOriginatingLevel const originating_level
):
	utility::pointer::ReferenceCount(),
	thread_request_originating_level_( originating_level )
{}

/// @brief Copy constructor.
RosettaThreadAssignmentInfo::RosettaThreadAssignmentInfo( RosettaThreadAssignmentInfo const & ) {}

/// @brief Destructor.
RosettaThreadAssignmentInfo::~RosettaThreadAssignmentInfo(){}

/// @brief Clone operation: copy this object and return an owning pointer to the copy.
RosettaThreadAssignmentInfoOP
RosettaThreadAssignmentInfo::clone() const {
	return utility::pointer::make_shared< RosettaThreadAssignmentInfo >( *this );
}

/// @brief Set the list of assigned threads (indices ranging from 1 to nthreads - 1 ).  These don't include the calling thread.
void
RosettaThreadAssignmentInfo::set_assigned_child_threads(
	utility::vector1< platform::Size > const & threadlist_in
) {
	//Obtain a write lock:
	utility::thread::WriteLockGuard lock( assigned_child_threads_mutex_ );
	assigned_child_threads_ = threadlist_in;
}

/// @brief Set the number of threads requested for the task.
void
RosettaThreadAssignmentInfo::set_requested_thread_count(
	platform::Size const setting
) {
	runtime_assert_string_msg( setting != 0, "Error in RosettaThreadAssignmentInfo::set_requested_thread_count(): The number of requested threads must be greater than zero." );
	//Obtain a write lock:
	utility::thread::WriteLockGuard lock( requested_thread_count_mutex_ );
	requested_thread_count_ = setting;
}

/// @brief Get the list of assigned threads (indices ranging from 1 to nthreads - 1 ).  These don't include the calling thread.
/// @details Deliberately passes result by copy rather than returning an instance, for thread-safety.
utility::vector1< platform::Size >
RosettaThreadAssignmentInfo::get_assigned_child_threads() const {
	//Obtain a read lock:
	utility::thread::ReadLockGuard lock( assigned_child_threads_mutex_ );
	return assigned_child_threads_;
}

/// @brief Get the number of additional threads that were assigned this task.  The total number of threads available for this
/// task are this number plus one, since the requesting thread will also run the task.
platform::Size
RosettaThreadAssignmentInfo::get_assigned_child_thread_count() const {
	//Obtain a read lock:
	utility::thread::ReadLockGuard lock( assigned_child_threads_mutex_ );
	return assigned_child_threads_.size();
}

/// @brief Get the total number of threads that were assigned this task.  This is one greater than the number of child threads
/// assigned this task, since the total includes the requesting thread, which also runs the task.
platform::Size
RosettaThreadAssignmentInfo::get_assigned_total_thread_count() const {
	//Obtain a read lock:
	utility::thread::ReadLockGuard lock( assigned_child_threads_mutex_ );
	return assigned_child_threads_.size() + 1;
}

/// @brief Get the number of threads that were originally requested for this task.  This is equal to or greater than the number
/// on which the task is actually running, since fewer threads might have been available.
platform::Size
RosettaThreadAssignmentInfo::get_requested_thread_count() const {
	//Obtain a read lock:
	utility::thread::ReadLockGuard lock( requested_thread_count_mutex_ );
	return requested_thread_count_;
}

/// @brief Get the level from which the thread request originally came.
RosettaThreadRequestOriginatingLevel
RosettaThreadAssignmentInfo::get_thread_request_originating_level() const {
	//Since this is set on object creation, there is no need for a lock to read this.
	return thread_request_originating_level_;
}

} //basic
} //thread_manager

#endif //MULTI_THREADED
