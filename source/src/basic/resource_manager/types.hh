// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/types.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_types_hh
#define INCLUDED_basic_resource_manager_types_hh

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <string>

namespace basic {
namespace resource_manager {

// A resource can be any type of rosetta-recongnized object
// After it is created, it will be const through its lifetime.
// You may very well desire that a resource should be constructed
// with very complex properties or data -- perhaps even data that
// depends on a Pose. That is fine: the Pose should also be a
// resource and should be requested from the ResourceManager
// during the construction of the object.
typedef utility::pointer::ReferenceCount Resource;
typedef utility::pointer::ReferenceCountCOP ResourceCOP;

} // namespace
} // namespace

#endif // include guard
