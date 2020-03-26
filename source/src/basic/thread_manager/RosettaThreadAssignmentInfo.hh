// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadAssignmentInfo.hh
/// @brief A class for storing information about the threads to which a function has been assigned by the RosettaThreadManager.
/// @details Accessor methods are fully threadsafe.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_basic_thread_manager_RosettaThreadAssignmentInfo_hh
#define INCLUDED_basic_thread_manager_RosettaThreadAssignmentInfo_hh

#ifdef MULTI_THREADED

#include <basic/thread_manager/RosettaThreadAssignmentInfo.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>
#include <utility/thread/ReadWriteMutex.hh>

// Platform headers
#include <platform/types.hh>

// C++ headers
#include <mutex>
#include <atomic>
#include <map>

namespace basic {
namespace thread_manager {

/// @brief A class for storing information about the threads to which a function has been assigned by the RosettaThreadManager.
class RosettaThreadAssignmentInfo : public utility::VirtualBase {

public:

	/// @brief Default constructor -- explicitly deleted.
	RosettaThreadAssignmentInfo() = delete;

	/// @brief Assignment constructor.
	RosettaThreadAssignmentInfo( RosettaThreadRequestOriginatingLevel const originating_level );

	/// @brief Copy constructor -- explicitly deleted.
	RosettaThreadAssignmentInfo(RosettaThreadAssignmentInfo const & src) = delete;

	/// @brief Destructor.
	virtual ~RosettaThreadAssignmentInfo();

	/// @brief Assignment operator -- explicitly deleted.
	RosettaThreadAssignmentInfo & operator = ( RosettaThreadAssignmentInfo const & ) = delete;

public: //Setters and getters

	/// @brief Set the list of assigned threads (indices ranging from 1 to nthreads - 1 ).  These don't include the calling thread.
	void set_assigned_child_threads( utility::vector1< platform::Size > const & threadlist_in );

	/// @brief Set the number of threads requested for the task.
	void set_requested_thread_count( platform::Size const setting );

	/// @brief Get the list of assigned threads (indices ranging from 1 to nthreads - 1 ).  These don't include the calling thread.
	/// @details Deliberately passes result by copy rather than returning an instance, for thread-safety.
	utility::vector1< platform::Size > get_assigned_child_threads() const;

	/// @brief Get the thread from which this task originated, included in the set of threads assigned to the task (but
	/// not in the assigned_child_threads_ list).
	inline platform::Size get_originating_thread() const { return originating_thread_; }

	/// @brief Get the number of additional threads that were assigned this task.  The total number of threads available for this
	/// task are this number plus one, since the requesting thread will also run the task.
	platform::Size get_assigned_child_thread_count() const;

	/// @brief Get the total number of threads that were assigned this task.  This is one greater than the number of child threads
	/// assigned this task, since the total includes the requesting thread, which also runs the task.
	platform::Size get_assigned_total_thread_count() const;

	/// @brief Get the number of threads that were originally requested for this task.  This is equal to or greater than the number
	/// on which the task is actually running, since fewer threads might have been available.
	platform::Size get_requested_thread_count() const;

	/// @brief Get the level from which the thread request originally came.
	RosettaThreadRequestOriginatingLevel get_thread_request_originating_level() const;

	/// @brief Checks the current thread's Rosetta index, and converts it into a one-based index in the vector of threads assigned to this
	/// task.  For example, if threads 4, 7, and 9 are assigned to this task, and thread 7 calls this function, it will return "2",
	/// since thread 7 is the 2nd of 3 threads assigned to this task.
	/// @details Assumes that calling thread is one that is assigned to this task!
	platform::Size get_this_thread_index_in_assigned_set() const;

private: //Setup

	/// @brief Set up the map of Rosetta thread index to one-based thread index for this task.
	/// @details For example, if this task has been assigned Rosetta threads 4, 7, and 9, then thread 4
	/// maps to assigned thread 1, thread 7 maps to assigned thread 2, and thread 9 maps to assigned thread
	/// 3.
	void set_up_rosetta_thread_index_to_assigned_set_index();

private: //Data

	/// @brief Mutex for accessing the list of threads to which a function has been assigned.
	mutable utility::thread::ReadWriteMutex assigned_child_threads_mutex_;

	/// @brief List of threads to which a function has been assigned.
	utility::vector1< platform::Size > assigned_child_threads_;

	/// @brief The thread from which this task originated, included in the set of threads assigned to the task (but
	/// not in the assigned_child_threads_ list).
	platform::Size originating_thread_ = 0;

	/// @brief Map of Rosetta thread index to one-based thread index for this task.
	/// @details For example, if this task has been assigned Rosetta threads 4, 7, and 9, then thread 4
	/// maps to assigned thread 1, thread 7 maps to assigned thread 2, and thread 9 maps to assigned thread
	/// 3.
	std::map< platform::Size, platform::Size > rosetta_thread_index_to_assigned_set_index_;

	/// @brief Number of threads requested.
	std::atomic< platform::Size > requested_thread_count_;

	/// @brief Rosetta level from which the therad request originally came.
	RosettaThreadRequestOriginatingLevel thread_request_originating_level_ = RosettaThreadRequestOriginatingLevel::UNKNOWN;

};


} //thread_manager
} //basic

#else //Not MULTI_THREADED

#include <basic/thread_manager/RosettaThreadAssignmentInfo.fwd.hh>
#include <utility/VirtualBase.hh>

namespace basic {
namespace thread_manager {

// Define empty class for non-multi-threaded build:
class RosettaThreadAssignmentInfo : public utility::VirtualBase {
public:
	RosettaThreadAssignmentInfo() = delete;
	RosettaThreadAssignmentInfo( RosettaThreadRequestOriginatingLevel const ) {}
	RosettaThreadAssignmentInfo( RosettaThreadAssignmentInfo const & ) = delete;
	~RosettaThreadAssignmentInfo() override = default;
};

}
}

#endif //MULTI_THREADED

#endif //INCLUDED_basic_thread_manager_RosettaThreadAssignmentInfo_hh
