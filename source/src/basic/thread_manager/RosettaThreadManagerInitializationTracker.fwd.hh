// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadManagerInitializationTracker.fwd.hh
/// @brief A singleton that tracks whether we have already launched threads or not.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_basic_thread_manager_RosettaThreadManagerInitializationTracker_fwd_hh
#define INCLUDED_basic_thread_manager_RosettaThreadManagerInitializationTracker_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace basic {
namespace thread_manager {

class RosettaThreadManagerInitializationTracker;

typedef utility::pointer::shared_ptr< RosettaThreadManagerInitializationTracker > RosettaThreadManagerInitializationTrackerOP;
typedef utility::pointer::shared_ptr< RosettaThreadManagerInitializationTracker const > RosettaThreadManagerInitializationTrackerCOP;

} //basic
} //thread_manager

#endif //INCLUDED_basic_thread_manager_RosettaThreadManagerInitializationTracker_fwd_hh
