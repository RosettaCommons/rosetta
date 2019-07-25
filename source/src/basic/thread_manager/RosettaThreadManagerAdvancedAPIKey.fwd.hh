// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadManagerAdvancedAPIKey.fwd.hh
/// @brief A class containing no private member data and only a constructor as a private member function, with friendship to only those classes that should be able
/// to access the advanced RosettaThreadManager API.  Since the advanced API requires an instance of a RosettaThreadManagerAdvancedAPIKey, this ensures that only
/// those classes can access the advanced API.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_basic_thread_manager_RosettaThreadManagerAdvancedAPIKey_fwd_hh
#define INCLUDED_basic_thread_manager_RosettaThreadManagerAdvancedAPIKey_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward declarations of classes that can instantiate the RosettaThreadManagerAdvancedAPIKey.

class RosettaThreadManagerTests;


// Forward
namespace basic {
namespace thread_manager {

class RosettaThreadManagerAdvancedAPIKey;

typedef utility::pointer::shared_ptr< RosettaThreadManagerAdvancedAPIKey > RosettaThreadManagerAdvancedAPIKeyOP;
typedef utility::pointer::shared_ptr< RosettaThreadManagerAdvancedAPIKey const > RosettaThreadManagerAdvancedAPIKeyCOP;

} //basic
} //thread_manager

#endif //INCLUDED_basic_thread_manager_RosettaThreadManagerAdvancedAPIKey_fwd_hh
