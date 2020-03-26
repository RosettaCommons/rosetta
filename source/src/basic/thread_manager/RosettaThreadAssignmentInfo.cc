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
#include <basic/thread_manager/RosettaThreadManager.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "basic.thread_manager.RosettaThreadAssignmentInfo" );


namespace basic {
namespace thread_manager {

/// @brief Default constructor.
RosettaThreadAssignmentInfo::RosettaThreadAssignmentInfo(
	RosettaThreadRequestOriginatingLevel const originating_level
):
	utility::VirtualBase(),
	requested_thread_count_(0),
	thread_request_originating_level_( originating_level )
{}

/// @brief Destructor.
RosettaThreadAssignmentInfo::~RosettaThreadAssignmentInfo(){}

/// @brief Set the list of assigned threads (indices ranging from 1 to nthreads - 1 ).  These don't include the calling thread.
void
RosettaThreadAssignmentInfo::set_assigned_child_threads(
	utility::vector1< platform::Size > const & threadlist_in
) {
	//Obtain a write lock:
	utility::thread::WriteLockGuard lock( assigned_child_threads_mutex_ );
	assigned_child_threads_ = threadlist_in;
	set_up_rosetta_thread_index_to_assigned_set_index();
}

/// @brief Set the number of threads requested for the task.
void
RosettaThreadAssignmentInfo::set_requested_thread_count(
	platform::Size const setting
) {
	runtime_assert_string_msg( setting != 0, "Error in RosettaThreadAssignmentInfo::set_requested_thread_count(): The number of requested threads must be greater than zero." );
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
	return requested_thread_count_;
}

/// @brief Get the level from which the thread request originally came.
RosettaThreadRequestOriginatingLevel
RosettaThreadAssignmentInfo::get_thread_request_originating_level() const {
	//Since this is set on object creation, there is no need for a lock to read this.
	return thread_request_originating_level_;
}

/// @brief Checks the current thread's Rosetta index, and converts it into a one-based index in the vector of threads assigned to this
/// task.  For example, if threads 4, 7, and 9 are assigned to this task, and thread 7 calls this function, it will return "2",
/// since thread 7 is the 2nd of 3 threads assigned to this task.
/// @details Assumes that calling thread is one that is assigned to this task!
platform::Size
RosettaThreadAssignmentInfo::get_this_thread_index_in_assigned_set() const {
	return rosetta_thread_index_to_assigned_set_index_.at( RosettaThreadManager::get_instance()->get_rosetta_thread_index() );
}

/// @brief Set up the map of Rosetta thread index to one-based thread index for this task.
/// @details For example, if this task has been assigned Rosetta threads 4, 7, and 9, then thread 4
/// maps to assigned thread 1, thread 7 maps to assigned thread 2, and thread 9 maps to assigned thread
/// 3.
void
RosettaThreadAssignmentInfo::set_up_rosetta_thread_index_to_assigned_set_index() {
	RosettaThreadManager * rtm( RosettaThreadManager::get_instance() );
	originating_thread_ = rtm->get_rosetta_thread_index();
	rosetta_thread_index_to_assigned_set_index_[originating_thread_] = 1;
	for ( platform::Size i(1), imax(assigned_child_threads_.size()); i<=imax; ++i ) {
		rosetta_thread_index_to_assigned_set_index_[assigned_child_threads_[i]] = i+1;
	}
}

} //basic
} //thread_manager

#endif //MULTI_THREADED
