// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadManagerInitializationTracker.cc
/// @brief A singleton that tracks whether we have already launched threads or not.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#include <basic/thread_manager/RosettaThreadManagerInitializationTracker.hh>


// Unit headers

// Project header

// Utility headers

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <string>

// Construct tracer.
static basic::Tracer TR( "basic.thread_manager.RosettaThreadManagerInitializationTracker" );

namespace basic {
namespace thread_manager {

// Public methods /////////////////////////////////////////////////////////////
// Static constant data access

// Private methods ////////////////////////////////////////////////////////////
// Empty constructor
RosettaThreadManagerInitializationTracker::RosettaThreadManagerInitializationTracker()
{}

/// @brief Determine whether thread launch has started.
bool
RosettaThreadManagerInitializationTracker::thread_manager_initialization_begun() const {
#ifdef MULTI_THREADED
	utility::thread::ReadLockGuard lock( initialization_mutex_ );
	return thread_manager_initialization_begun_;
#else
	return false; //Always false in non-multi-threaded build.
#endif
}

/// @brief Store the fact that thread lauch has started.
void
RosettaThreadManagerInitializationTracker::mark_thread_manager_initialization_as_begun() {
#ifdef MULTI_THREADED
	utility::thread::WriteLockGuard lock( initialization_mutex_ );
#endif
	thread_manager_initialization_begun_ = true;
}

/// @brief Determine whether threads have been launched.
bool
RosettaThreadManagerInitializationTracker::thread_manager_was_initialized() const {
#ifdef MULTI_THREADED
	utility::thread::ReadLockGuard lock( initialization_mutex_ );
	return thread_manager_was_initialized_;
#else
	return false; //Always false in non-multi-threaded build.
#endif
}

/// @brief Store the fact that threads have been launched.
void
RosettaThreadManagerInitializationTracker::mark_thread_manager_as_initialized() {
#ifdef MULTI_THREADED
	utility::thread::WriteLockGuard lock( initialization_mutex_ );
#endif
	thread_manager_was_initialized_ = true;
}

} //basic
} //thread_manager
