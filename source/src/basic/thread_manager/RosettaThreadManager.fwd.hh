// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadManager.fwd.hh
/// @brief A manager that maintains a threadpool and handles requests for threads for multithreaded
/// execution of functions.  This allows multithreading at many different levels in the Rosetta
/// library hierarchy, from job-level parallel execution down to parallel computation of a score,
/// gradient vector, or interaction graph.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef MULTI_THREADED

#ifndef INCLUDED_basic_thread_manager_RosettaThreadManager_fwd_hh
#define INCLUDED_basic_thread_manager_RosettaThreadManager_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

// C++ headers
#include <functional>

// Forward
namespace basic {
namespace thread_manager {

class RosettaThreadManager;

typedef utility::pointer::shared_ptr< RosettaThreadManager > RosettaThreadManagerOP;
typedef utility::pointer::shared_ptr< RosettaThreadManager const > RosettaThreadManagerCOP;
typedef utility::pointer::weak_ptr< RosettaThreadManager > RosettaThreadManagerAP;
typedef utility::pointer::weak_ptr< RosettaThreadManager const > RosettaThreadManagerCAP;

typedef std::function< void () > RosettaThreadFunction;
typedef utility::pointer::shared_ptr< std::function< void () > > RosettaThreadFunctionOP;

} //thread_manager
} //basic

#endif //INCLUDED_basic_thread_manager_RosettaThreadManager_fwd_hh

#endif //MULTI_THREADED
