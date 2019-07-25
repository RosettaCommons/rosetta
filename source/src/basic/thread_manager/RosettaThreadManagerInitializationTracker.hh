// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadManagerInitializationTracker
/// @brief A singleton that tracks whether we have already launched threads or not.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_basic_thread_manager_RosettaThreadManagerInitializationTracker_hh
#define INCLUDED_basic_thread_manager_RosettaThreadManagerInitializationTracker_hh

// Unit headers
#include <basic/thread_manager/RosettaThreadManagerInitializationTracker.fwd.hh>

// Utility header
#include <utility/SingletonBase.hh>

#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#endif

// C++ header

namespace basic {
namespace thread_manager {

/// @brief A singleton that tracks whether we have already launched threads or not.
class RosettaThreadManagerInitializationTracker : public utility::SingletonBase< RosettaThreadManagerInitializationTracker > {
	friend class utility::SingletonBase< RosettaThreadManagerInitializationTracker >;

private:  // Private methods //////////////////////////////////////////////////
	// Empty constructor
	RosettaThreadManagerInitializationTracker();
	RosettaThreadManagerInitializationTracker(RosettaThreadManagerInitializationTracker const & ) = delete;
	RosettaThreadManagerInitializationTracker operator=(RosettaThreadManagerInitializationTracker const & ) = delete;

public:

	/// @brief Determine whether thread launch has started.
	bool thread_manager_initialization_begun() const;

	/// @brief Store the fact that thread lauch has started.
	void mark_thread_manager_initialization_as_begun();

	/// @brief Determine whether threads have been launched.
	bool thread_manager_was_initialized() const;

	/// @brief Store the fact that threads have been launched.
	void mark_thread_manager_as_initialized();

private:  // Private data /////////////////////////////////////////////////////

	bool thread_manager_initialization_begun_ = false;

	bool thread_manager_was_initialized_ = false;

#ifdef MULTI_THREADED
	mutable utility::thread::ReadWriteMutex initialization_mutex_;
#endif

};

} //basic
} //thread_manager

#endif //INCLUDED_basic/thread_manager_RosettaThreadManagerInitializationTracker_fwd_hh



