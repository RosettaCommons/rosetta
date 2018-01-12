// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/ResourceManagerCreator.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_manager_ResourceManagerCreator_FWD_HH
#define INCLUDED_basic_resource_manager_ResourceManagerCreator_FWD_HH

// package headers
#include <basic/resource_manager/ResourceManager.fwd.hh>

//C++ headers
#include <string>

namespace basic {
namespace resource_manager {

class ResourceManagerCreator;
typedef utility::pointer::shared_ptr< ResourceManagerCreator > ResourceManagerCreatorOP;
typedef utility::pointer::shared_ptr< ResourceManagerCreator const > ResourceManagerCreatorCOP;

} // namespace resource_manager
} // namespace basic

#endif //INCLUDED_basic_resource_manager_ResourceManager_HH
