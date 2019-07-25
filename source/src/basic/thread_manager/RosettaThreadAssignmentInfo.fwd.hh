// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadAssignmentInfo.fwd.hh
/// @brief A class for storing information about the threads to which a function has been assigned by the RosettaThreadManager.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_basic_thread_manager_RosettaThreadAssignmentInfo_fwd_hh
#define INCLUDED_basic_thread_manager_RosettaThreadAssignmentInfo_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace basic {
namespace thread_manager {

enum class RosettaThreadRequestOriginatingLevel {
	UNKNOWN = 1, //Keep first
	CORE_GENERIC,
	CORE_PACK,
	PROTOCOLS_GENERIC,
	PROTOCOLS_JOB_DISTRIBUTOR,
	APPLICATIONS_OR_APPLICATION_PROTOCOLS,
	INVALID, //Keep second-to-last
	END_OF_LIST = INVALID //Keep last
};

class RosettaThreadAssignmentInfo;

typedef utility::pointer::shared_ptr< RosettaThreadAssignmentInfo > RosettaThreadAssignmentInfoOP;
typedef utility::pointer::shared_ptr< RosettaThreadAssignmentInfo const > RosettaThreadAssignmentInfoCOP;

} //thread_manager
} //basic

#endif //INCLUDED_basic_thread_manager_RosettaThreadAssignmentInfo_fwd_hh
