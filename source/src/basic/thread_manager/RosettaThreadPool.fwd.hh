// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadPool.fwd.hh
/// @brief A container for a pool of threads.  Threads idle continuously until loaded with a
/// function to execute.  A subset of threads then execute the function synchronously and then
/// return to the idle state.
/// @note The basic verison of this class assigns functions to threads on a first-come, first-served
/// basis.  If I request 8 threads and 8 threads are idle, I get all 8 threads.  If the next request
/// also asks for 8 threads while I'm still using mine, it only gets the thread that made the request.
/// Derived classes can implement different behaviour for assigning work to threads.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef MULTI_THREADED

#ifndef INCLUDED_basic_thread_manager_RosettaThreadPool_fwd_hh
#define INCLUDED_basic_thread_manager_RosettaThreadPool_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>


// Forward
namespace basic {
namespace thread_manager {

class RosettaThreadPool;

typedef utility::pointer::shared_ptr< RosettaThreadPool > RosettaThreadPoolOP;
typedef utility::pointer::shared_ptr< RosettaThreadPool const > RosettaThreadPoolCOP;
typedef utility::pointer::weak_ptr< RosettaThreadPool > RosettaThreadPoolAP;
typedef utility::pointer::weak_ptr< RosettaThreadPool const > RosettaThreadPoolCAP;

} //thread_manager
} //basic

#endif //INCLUDED_basic_thread_manager_RosettaThreadPool_fwd_hh

#endif //MULTI_THREADED
