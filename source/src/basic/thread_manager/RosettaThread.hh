// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThread.hh
/// @brief A container for a thread in a RosettaThreadPool.  The thread idles continuously until
/// loaded with a function to execute.  It then executes the function and returns to the idle state.
/// @note This object is held by the RosettaThreadPool, which is intended to be held by the global
/// static singleton RosettaThreadManager.  As such, it can potentially be accessed from any thread.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef MULTI_THREADED

#ifndef INCLUDED_basic_thread_manager_RosettaThread_hh
#define INCLUDED_basic_thread_manager_RosettaThread_hh

#include <basic/thread_manager/RosettaThread.fwd.hh>
#include <basic/thread_manager/RosettaThreadPool.fwd.hh>
#include <basic/thread_manager/RosettaThreadManager.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// Platform headers
#include <platform/types.hh>

// STL headers
#include <functional>
#include <mutex>
#include <thread>
#include <atomic>
#include <condition_variable>

namespace basic {
namespace thread_manager {

/// @brief A "key" class that is used to ensure that ONLY the RosettaThreadPool can instantiate the RosettaThread class.
/// @note This class declares the RosettaThreadPool as a friend.  Since this class has NO private data and NO private members (aside
/// from the constructor), this does not jeopardize const correctness in any way.
class RosettaThreadInstantiationKey {

	friend class RosettaThreadPool;

private:

	/// @brief Constructor.  Private, to ensure that only the RosettaThreadManager can instantiate this.
	RosettaThreadInstantiationKey() = default;

	/// @brief Copy constructor.  Private, to ensure that only the RosettaThreadManager can copy this.
	RosettaThreadInstantiationKey( RosettaThreadInstantiationKey const & ) = default;

public:

	/// @brief Destructor.
	~RosettaThreadInstantiationKey() = default;

};

/// @brief A container for a thread in a RosettaThreadPool.  The thread idles continuously until loaded with a function to execute.
/// It then executes the function and returns to the idle state.
class RosettaThread : public utility::pointer::ReferenceCount {

public:

	/// @brief Default constructor -- explicitly deleted.
	RosettaThread() = delete;

	/// @brief RosettaThread constructor.  The thread index *must* be passed.
	/// @note A RosettaThreadInstantiationKey must also be passed to prove that this thread is being launched
	/// ONLY in a permitted context (the context of the RosettaThreadPool).
	RosettaThread( platform::Size const thread_index, RosettaThreadInstantiationKey const & key );

	/// @brief Copy constructor -- explicitly deleted.
	RosettaThread(RosettaThread const &) = delete;

	/// @brief Destructor.  Calls force_thread_termination().
	virtual ~RosettaThread();

public:

	/// @brief Set the function to run on this thread.
	/// @details Requires a function that returns void (bundled with its arguments with std::bind).
	/// @returns False for failure (thread was occupied), true for success (function was assigned to thread).
	/// @note The job_completion_* objects are used to signal job completion to the RosettaThreadPool.
	bool set_function( RosettaThreadFunctionOP function_to_execute, std::mutex & job_completion_mutex, std::condition_variable & job_completion_condition_variable, platform::Size & job_completion_count );

	/// @brief Set the idle state of this thread.  If set to true and a function is already running, it continues
	/// to run, but no new function will run until this is set to false.
	void set_forced_idle( bool const setting );

	/// @brief Is this thread idle?
	/// @details This function locks the thread mutex, and should not be called from a locked context.
	bool is_idle() const;

	/// @brief Indicate that this thread is available to be assigned new work.
	inline void set_available_for_new_work() { is_available_for_new_work_ = true; /*Atomic write operation.*/ }

	/// @brief Check whether this thread is available to be assigned new work.
	inline bool is_available_for_new_work() const { return is_available_for_new_work_.load(); /*Atomic read operation.*/ }

	/// @brief Get the system ID of this thread.
	inline std::thread::id get_system_thread_id() const { return this_thread_.get_id(); }

	/// @brief Get the thread index, in Rosetta numbering.  Thread 0 is the master; subsequent threads are numbered from zero.
	inline platform::Size get_thread_index() const { return thread_index_; }

private:

	/// @brief The function that this thread performs, which idles until the function object is populated.
	void thread_function();

	/// @brief Terminates the std::thread object held by this object.
	/// @details Warning!  This should only be called by this object's destructor.  Renders the object unusable afterwards.
	void force_thread_termination();

	/// @brief Checks whether there's a function to execute.  Returns true if there is NOT.
	/// @details This function does NOT lock the thread mutex.  It should only be called in a context in which
	/// the thread mutex is already locked.
	inline bool
	is_idle_already_locked() const {
		return function_to_execute_ == nullptr;
	}

private:

	/// @brief Index of this thread, in Rosetta numbering.  The master thread is thread zero, so this must be greater than zero.
	platform::Size thread_index_ = 0;

	/// @brief A mutex for locking the thread for loading a function to run.
	mutable std::mutex this_thread_function_mutex_;

	/// @brief A condition variable, used to wake this thread up when there is work to be done.
	std::condition_variable cv_;

	/// @brief If true, the thread is prevented from running anything.  False by default.  Can be temporarily set to true
	/// by external code.
	std::atomic_bool hold_thread_;

	/// @brief If true, then this thread can have new work assigned to it.
	std::atomic_bool is_available_for_new_work_;

	/// @brief The function to execute.
	RosettaThreadFunctionOP function_to_execute_ = nullptr;

	/// @brief The thread managed by this object.
	std::thread this_thread_;

	/// @brief Boolean used to force termination of a thread.
	std::atomic_bool terminate_now_;

	/// @brief A pointer to a mutex used to signal job completion.
	/// @note This must be a pointer and not a reference since it will point to different objects at different times.  The original
	/// object was not created by smart pointer due to the undue overhead associated with that.  It is held by the RosettaThreadPool,
	/// and is guaranteed to exist for the duration of the time that this pointer is non-null.
	std::mutex * job_completion_mutex_ = nullptr;

	/// @brief A pointer to a condition variable used to signal job completion.
	/// @note This must be a pointer and not a reference since it will point to different objects at different times.  The original
	/// object was not created by smart pointer due to the undue overhead associated with that.  It is held by the RosettaThreadPool,
	/// and is guaranteed to exist for the duration of the time that this pointer is non-null.
	std::condition_variable * job_completion_condition_variable_ = nullptr;

	/// @brief A pointer to a counter used to count the number of threads that were assigned a job that have
	/// completed the job.  Incremented by one by this thread on job completion.
	/// @note This must be a pointer and not a reference since it will point to different objects at different times.  The original
	/// object was not created by smart pointer due to the undue overhead associated with that.  It is held by the RosettaThreadPool,
	/// and is guaranteed to exist for the duration of the time that this pointer is non-null.
	platform::Size * job_completion_count_ = nullptr;

};


} //thread_manager
} //basic

#endif //INCLUDED_basic_thread_manager_RosettaThread_hh
#endif //MULTI_THREADED
