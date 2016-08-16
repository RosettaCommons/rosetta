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

// A resource can be any type of rosetta recongnized object
typedef utility::pointer::ReferenceCount Resource;
typedef utility::pointer::ReferenceCountOP ResourceOP;

// This identifies a resource, eg native, fragments:3mer, native:electron_density
// It may make sense to use the utility::key framework for this
typedef std::string ResourceDescription;

// These identify specific instances of a resource, 1j2n
typedef std::string ResourceTag;

// This identifies which instance of a ResourceDescription should be returned to a protocol
typedef std::string JobTag;

// Identify a resource locator with the ResourceLocatorFactory
typedef std::string LocatorType;
typedef std::string LocatorTag;
typedef std::string LocatorID;


// Identify a resource loader with the ResourceLoaderFactory
typedef std::string LoaderType;

typedef std::string ResourceOptionsTag;


} // namespace
} // namespace

#endif // include guard
