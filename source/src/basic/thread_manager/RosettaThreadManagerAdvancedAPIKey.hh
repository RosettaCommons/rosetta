// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadManagerAdvancedAPIKey.hh
/// @brief A class containing no private member data and only a constructor as a private member function, with friendship to only those classes that should be able
/// to access the advanced RosettaThreadManager API.  Since the advanced API requires an instance of a RosettaThreadManagerAdvancedAPIKey, this ensures that only
/// those classes can access the advanced API.'
///
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
///
/// @details If you use the RosettaThreadManager::run_function_in_threads() advanced API, you must modify this class to add your calling context to the friend list
/// for this class.  If you do not do this, you will not be able to create an instance of this class to pass to that function.  Whitelisting your calling context
/// will light up the central_class_modification integration test, and this means that you will be forced to justify why you need the advanced API before you can
/// merge your code.  It is STRONGLY recommended that you try to use the basic API (the RosettaThreadManager::do_work_vector_in_threads() function), which is much
/// safer, if you possibly can.  If you want to discuss the best way to multi-thread your code, please talk to Andrew Leaver-Fay or Vikram K. Mulligan BEFORE
/// starting.


#ifndef INCLUDED_basic_thread_manager_RosettaThreadManagerAdvancedAPIKey_hh
#define INCLUDED_basic_thread_manager_RosettaThreadManagerAdvancedAPIKey_hh

#include <basic/thread_manager/RosettaThreadManagerAdvancedAPIKey.fwd.hh>
#include <basic/thread_manager/RosettaThreadManager.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace basic {
namespace thread_manager {

/// @brief A class containing no private member data and only a constructor as a private member function, with friendship to only those classes that should be able
/// to access the advanced RosettaThreadManager API.  Since the advanced API requires an instance of a RosettaThreadManagerAdvancedAPIKey, this ensures that only
/// those classes can access the advanced API.
///
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
///
/// @details If you use the RosettaThreadManager::run_function_in_threads() advanced API, you must modify this class to add your calling context to the friend list
/// for this class.  If you do not do this, you will not be able to create an instance of this class to pass to that function.  Whitelisting your calling context
/// will light up the central_class_modification integration test, and this means that you will be forced to justify why you need the advanced API before you can
/// merge your code.  It is STRONGLY recommended that you try to use the basic API (the RosettaThreadManager::do_work_vector_in_threads() function), which is much
/// safer, if you possibly can.  If you want to discuss the best way to multi-thread your code, please talk to Andrew Leaver-Fay or Vikram K. Mulligan BEFORE
/// starting.
class RosettaThreadManagerAdvancedAPIKey {

	// Whitelist classes here.  Note that whitelisting a class to access the advanced RosettaThreadManager API
	// requires a justification!
	friend class RosettaThreadManager; // The basic API calls the advanced API.
	friend class ::RosettaThreadManagerTests; //To allow testing.

private:

	/// @brief Private, explicitly-declared constructor ensures that we can't instantiate except in whitelisted contexts.
	RosettaThreadManagerAdvancedAPIKey() {};

	/// @brief Prohibit copying.
	RosettaThreadManagerAdvancedAPIKey( RosettaThreadManagerAdvancedAPIKey const & ) = delete;
public:
	~RosettaThreadManagerAdvancedAPIKey() = default;
};

} //basic
} //thread_manager

#endif //INCLUDED_basic_thread_manager_RosettaThreadManagerAdvancedAPIKey_hh
