// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadPool.hh
/// @brief A container for a pool of threads.  Threads idle continuously until loaded with a
/// function to execute.  A subset of threads then execute the function synchronously and then
/// return to the idle state.
/// @note The basic verison of this class assigns functions to threads on a first-come, first-served
/// basis.  If I request 8 threads and 8 threads are idle, I get all 8 threads.  If the next request
/// also asks for 8 threads while I'm still using mine, it only gets the thread that made the request.
/// Derived classes can implement different behaviour for assigning work to threads.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_basic_thread_manager_RosettaThreadPool_hh
#define INCLUDED_basic_thread_manager_RosettaThreadPool_hh

#ifdef MULTI_THREADED

#include <basic/thread_manager/RosettaThreadPool.fwd.hh>
#include <basic/thread_manager/RosettaThread.fwd.hh>
#include <basic/thread_manager/RosettaThreadAssignmentInfo.fwd.hh>
#include <basic/thread_manager/RosettaThreadManager.fwd.hh>
#include <basic/thread_manager/RosettaThreadAllocation.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

// STL headers
#include <functional>
#include <mutex>
#include <thread>
#include <map>

// Platform headers
#include <platform/types.hh>

#include <utility/vector1.hh>

namespace basic {
namespace thread_manager {

/// @brief A "key" class that is used to ensure that ONLY the RosettaThreadManager can instantiate the RosettaThreadPool class.
/// @note This class declares the RosettaThreadManager as a friend.  Since this class has NO private data and NO private members (aside
/// from the constructor), this does not jeopardize const correctness in any way.
class RosettaThreadPoolInstantiationKey {

	friend class RosettaThreadManager;

private:

	/// @brief Constructor.  Private, to ensure that only the RosettaThreadManager can instantiate this.
	RosettaThreadPoolInstantiationKey() = default;

	/// @brief Copy constructor.  Private, to ensure that only the RosettaThreadManager can copy this.
	RosettaThreadPoolInstantiationKey( RosettaThreadPoolInstantiationKey const & ) = default;

public:

	/// @brief Destructor.
	~RosettaThreadPoolInstantiationKey() = default;

};

/// @brief A container for a pool of threads.  Threads idle continuously until loaded with a function to execute.  A subset of
/// threads then execute the function synchronously and then return to the idle state.
class RosettaThreadPool : public utility::VirtualBase {

public:

	/// @brief Default constructor -- explicitly deleted.
	RosettaThreadPool() = delete;

	/// @brief Constructor requiring the number of threads to launch.
	/// @note The pool holds total_threads-1 threads.  A RosettaThreadPoolInstantiationKey must be passed to prove that
	/// this is being instantiated in a context that is allowed to instantiate the RosettaThreadPool (i.e. the context of
	/// the RosettaThreadManager).
	RosettaThreadPool( platform::Size const total_threads, RosettaThreadPoolInstantiationKey const & key );

	/// @brief Copy constructor -- explicitly deleted.
	RosettaThreadPool(RosettaThreadPool const &) = delete;

	/// @brief Destructor.
	/// @details Calls force_terminate_threads();
	virtual ~RosettaThreadPool();

public: //Functions

	/// @brief Given a function that was bundled with its arguments with std::bind, run it in many threads.
	///
	/// @details The function is assigned to as many threads as are available, including the thread from which the request
	/// originates.  It is guaranteed to run in 1 <= actual_thread_count <= requested_thread_count threads.  After assigning
	/// the function to up to (requsted_thread_count - 1) other threads, the function executes in the current thread, then
	/// the current thread blocks until the assigned threads report that they are idle.
	///
	/// The RosettaThreadManager must provide a strong guarantee that the RosettaThreadFunction object exists and will continue
	/// to exist until execution terminates.  Since the RosettaThreadManager is the only class that can instantiate or access
	/// the RosettaThreadPool, this is a safe ownership arrangement.
	///
	/// @note A RosettaThreadAssignmentInfo object should be passed in.  It will be populated with
	/// the number of threads requested, the number actually assigned, the indices of the assigned threads, and a map of
	/// system thread ID to Rosetta thread index.  The same owning pointer may optionally be provided to the function to
	/// execute by the calling function if the function to execute requires access to this information.
	virtual
	void
	run_function_in_threads(
		RosettaThreadFunction * function_to_execute,
		platform::Size const requested_thread_count,
		RosettaThreadAssignmentInfo & thread_assignment
	);


	//PREALLOCATION SECTION

	RosettaThreadAllocation
	preallocate_threads(
		platform::Size const requested_thread_count,
		RosettaThreadAssignmentInfo & thread_assignment
	);

	void
	release_threads(
		RosettaThreadAllocation & allocation
	);

	void
	run_function_in_threads(
		RosettaThreadFunction * function_to_execute,
		RosettaThreadAllocation & allocation
	);


	/// @brief Force all threads to terminate.  Called by the destructor of the RosettaThreadManager class, and
	/// by the destructor of this class.
	/// @details At the end of this operation, the threads_ object is empty.
	void
	force_terminate_threads();

	/// @brief Construct a map of (system thread id)->(Rosetta thread index), and return the object.
	/// @details Intended to be called only by the RosettaThreadManager.  If you want the current thread ID, use
	/// basic::thread_manager::RosettaThreadManager::get_instance()->get_rosetta_thread_index().
	std::map< std::thread::id, platform::Size > get_thread_id_to_rosetta_thread_index_map() const;

private: //Data

	/// @brief The pool of threads, created on object initialization, destroyed on
	/// destruction, and kept running in between.
	utility::vector1< RosettaThreadOP > threads_;

	/// @brief A mutex for locking the list of threads that are working.
	mutable std::mutex working_thread_list_mutex_;

};

} //thread_manager
} //basic

#endif //MULTI_THREADED
#endif //INCLUDED_basic_thread_manager_RosettaThreadPool_hh
