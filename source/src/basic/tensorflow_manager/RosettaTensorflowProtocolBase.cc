// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  basic/tensorflow_manager/RosettaTensorflowProtocolBase.cc
/// @brief A pure virtual base class for storing Tensorflow sessions and the code for executing them and returning a result.
/// Derived classes will be protocol-specific, and will accept a RosettaTensorflowInput and produce a RosettaTensorflowOutput.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Unit header or inline function header
#include <basic/tensorflow_manager/RosettaTensorflowProtocolBase.hh>

// NOTE: This file should have NO dependencies other than its header.


namespace basic {
namespace tensorflow_manager {


// Defined to prevent pure virtual destructor error at run time.
RosettaTensorflowProtocolBase::~RosettaTensorflowProtocolBase(){}


} //tensorflow_manager
} //basic


