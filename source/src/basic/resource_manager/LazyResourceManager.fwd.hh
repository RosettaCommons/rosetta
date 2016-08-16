// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/LazyResourceManager.fwd.hh
/// @brief  Forward Header for base class for LazyResourceManager
/// @author Matthew O'Meara

#ifndef INCLUDED_basic_resource_manager_LazyResourceManager_fwd_hh
#define INCLUDED_basic_resource_manager_LazyResourceManager_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace basic {
namespace resource_manager {

class LazyResourceManager;
typedef utility::pointer::shared_ptr< LazyResourceManager > LazyResourceManagerOP;
typedef utility::pointer::shared_ptr< LazyResourceManager const > LazyResourceManagerCOP;

} //namespace
} //namespace

#endif
