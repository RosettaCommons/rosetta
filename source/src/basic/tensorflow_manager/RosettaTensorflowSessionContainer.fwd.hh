// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/tensorflow_manager/RosettaTensorflowSessionContainer.fwd.hh
/// @brief A container for Rosetta Tensorflow sessions, allowing sessions to be loaded once and stored in the global Tensorflow session manager.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_basic_tensorflow_manager_RosettaTensorflowSessionContainer_fwd_hh
#define INCLUDED_basic_tensorflow_manager_RosettaTensorflowSessionContainer_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace basic {
namespace tensorflow_manager {

class RosettaTensorflowSessionContainer;

using RosettaTensorflowSessionContainerOP = utility::pointer::shared_ptr< RosettaTensorflowSessionContainer >;
using RosettaTensorflowSessionContainerCOP = utility::pointer::shared_ptr< RosettaTensorflowSessionContainer const >;

} //tensorflow_manager
} //basic

#endif //INCLUDED_basic_tensorflow_manager_RosettaTensorflowSessionContainer_fwd_hh
