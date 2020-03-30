// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/tensorflow_manager/RosettaTensorflowProtocolBase.fwd.hh
/// @brief A pure virtual base class for storing Tensorflow sessions and the code for executing them and returning a result.
/// Derived classes will be protocol-specific, and will accept a RosettaTensorflowInput and produce a RosettaTensorflowOutput.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifndef INCLUDED_basic_tensorflow_manager_RosettaTensorflowProtocolBase_fwd_hh
#define INCLUDED_basic_tensorflow_manager_RosettaTensorflowProtocolBase_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace basic {
namespace tensorflow_manager {

class RosettaTensorflowProtocolBase;

using RosettaTensorflowProtocolBaseOP = utility::pointer::shared_ptr< RosettaTensorflowProtocolBase >;
using RosettaTensorflowProtocolBaseCOP = utility::pointer::shared_ptr< RosettaTensorflowProtocolBase const >;

} //tensorflow_manager
} //basic

#endif //INCLUDED_basic_tensorflow_manager_RosettaTensorflowProtocolBase_fwd_hh
