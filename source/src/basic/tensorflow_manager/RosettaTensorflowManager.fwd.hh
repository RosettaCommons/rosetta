// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/tensorflow_manager/RosettaTensorflowManager.fwd.hh
/// @brief A manager class for loading Tensorflow sessions and controlling mapping to the CPU or GPU.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_basic_tensorflow_manager_RosettaTensorflowManager_fwd_hh
#define INCLUDED_basic_tensorflow_manager_RosettaTensorflowManager_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// STL headers
#include <chrono>

// Forward
namespace basic {
namespace tensorflow_manager {

#if defined(MAC) || defined(__APPLE__)  ||  defined(__OSX__)
typedef std::chrono::system_clock ROSETTA_TENSORFLOW_CLOCK;
#else
typedef std::chrono::high_resolution_clock ROSETTA_TENSORFLOW_CLOCK;
#endif

class RosettaTensorflowManager;

using RosettaTensorflowManagerOP = utility::pointer::shared_ptr< RosettaTensorflowManager >;
using RosettaTensorflowManagerCOP = utility::pointer::shared_ptr< RosettaTensorflowManager const >;

} //tensorflow_manager
} //basic

#endif //INCLUDED_basic_tensorflow_manager_RosettaTensorflowManager_fwd_hh
