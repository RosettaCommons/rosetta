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
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/thread/ReadWriteMutex.hh>

// Platform headers
#include <platform/types.hh>

// C++ headers
#include <mutex>

namespace basic {
namespace thread_manager {

/// @brief A class for storing information about the threads to which a function has been assigned by the RosettaThreadManager.
class RosettaThreadAssignmentInfo : public utility::pointer::ReferenceCount {

public:

	/// @brief Default constructor -- explicitly deleted.
	RosettaThreadAssignmentInfo() = delete;

	/// @brief Assignment constructor.
	RosettaThreadAssignmentInfo( RosettaThreadRequestOriginatingLevel const originating_level );

	/// @brief Copy constructor.
	RosettaThreadAssignmentInfo(RosettaThreadAssignmentInfo const & src);

	/// @brief Destructor.
	virtual ~RosettaThreadAssignmentInfo();

	/// @brief Clone operation: copy this object and return an owning pointer to the copy.
	RosettaThreadAssignmentInfoOP clone() const;

public: //Setters and getters

	/// @brief Set the list of assigned threads (indices ranging from 1 to nthreads - 1 ).  These don't include the calling thread.
	void set_assigned_child_threads( utility::vector1< platform::Size > const & threadlist_in );

	/// @brief Set the number of threads requested for the task.
	void set_requested_thread_count( platform::Size const setting );

	/// @brief Get the list of assigned threads (indices ranging from 1 to nthreads - 1 ).  These don't include the calling thread.
	/// @details Deliberately passes result by copy rather than returning an instance, for thread-safety.
	utility::vector1< platform::Size > get_assigned_child_threads() const;

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

private: //Data

	/// @brief Mutex for accessing the list of threads to which a function has been assigned.
	mutable utility::thread::ReadWriteMutex assigned_child_threads_mutex_;

	/// @brief List of threads to which a function has been assigned.
	utility::vector1< platform::Size > assigned_child_threads_;

	/// @brief Mutex for accessing the count of requested threads.
	mutable utility::thread::ReadWriteMutex requested_thread_count_mutex_;

	/// @brief Number of threads requested.
	platform::Size requested_thread_count_ = 0;

	RosettaThreadRequestOriginatingLevel thread_request_originating_level_ = RosettaThreadRequestOriginatingLevel::UNKNOWN;

};


} //thread_manager
} //basic

#else //Not MULTI_THREADED

#include <basic/thread_manager/RosettaThreadAssignmentInfo.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace basic {
namespace thread_manager {

// Define empty class for non-multi-threaded build:
class RosettaThreadAssignmentInfo : public utility::pointer::ReferenceCount {
public:
	RosettaThreadAssignmentInfo() = delete;
	RosettaThreadAssignmentInfo( RosettaThreadRequestOriginatingLevel const ) {}
	RosettaThreadAssignmentInfo( RosettaThreadAssignmentInfo const & ) = default;
	~RosettaThreadAssignmentInfo() = default;
};

}
}

#endif //MULTI_THREADED

#endif //INCLUDED_basic_thread_manager_RosettaThreadAssignmentInfo_hh
