// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/locator/FileSystemResourceLocator.fwd.hh
/// @brief  forward header for FileSystemResourceLocator

#ifndef INCLUDED_basic_resource_manager_locator_FileSystemResourceLocator_fwd_hh
#define INCLUDED_basic_resource_manager_locator_FileSystemResourceLocator_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace basic {
namespace resource_manager {
namespace locator {

class FileStream;
typedef utility::pointer::shared_ptr< FileStream > FileStreamOP;
typedef utility::pointer::shared_ptr< FileStream const > FileStreamCOP;

class FileSystemResourceLocator;
typedef utility::pointer::shared_ptr< FileSystemResourceLocator > FileSystemResourceLocatorOP;
typedef utility::pointer::shared_ptr< FileSystemResourceLocator const > FileSystemResourceLocatorCOP;

}// namespace
}// namespace
}// namespace

#endif // include guard
