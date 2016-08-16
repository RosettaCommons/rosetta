// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/locator/NullResourceLocator.fwd.hh
/// @brief  forward header for NullResourceLocator

#ifndef INCLUDED_basic_resource_manager_locator_NullResourceLocator_fwd_hh
#define INCLUDED_basic_resource_manager_locator_NullResourceLocator_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace basic {
namespace resource_manager {
namespace locator {

class NullStream;
typedef utility::pointer::shared_ptr< NullStream > NullStreamOP;
typedef utility::pointer::shared_ptr< NullStream const > NullStreamCOP;

class NullResourceLocator;
typedef utility::pointer::shared_ptr< NullResourceLocator > NullResourceLocatorOP;
typedef utility::pointer::shared_ptr< NullResourceLocator const > NullResourceLocatorCOP;

}// namespace
}// namespace
}// namespace

#endif // include guard
