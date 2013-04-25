// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/FallbackConfiguration.hh
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_basic_resource_manager_fallback_configuration_HH
#define INCLUDED_basic_resource_manager_fallback_configuration_HH

// Unit Headers
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>

// Package headers
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>


//utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace basic {
namespace resource_manager {

/// @brief This class describes how a resource should be created
/// if the ResourceManager is not being used, e.g., if resources
/// have been specified on the command line.  A protocol
/// will still query the ResourceManager for its desired resource
/// (by name) and if no resource has been provided for that
/// description, then the ResourceManager will ask the
/// FallbackConfiguration for that Resource.
class FallbackConfiguration : public utility::pointer::ReferenceCount {

public:
	FallbackConfiguration() {}

	virtual ~FallbackConfiguration();

	/// @brief Has a fallback been provided for the given resource description?
	/// It is possible that the fallback configuration would be unable to deliver
	/// a requested resource.
	virtual
	bool
	fallback_specified( ResourceDescription const & desc ) const = 0;

	/// @brief What kind of resource loader should be used to instantiate the
	/// desired resource?
	virtual
	LoaderType
	get_resource_loader( ResourceDescription const & desc ) const = 0;

	/// @brief What is the locator tag that the resource locator should use to
	/// find the data used to construct the desired resource?
	virtual
	LocatorID
	get_locator_id( ResourceDescription const & desc ) const = 0;

	/// @brief Return a pointer to the resource options object that the
	/// ResourceLoader will use to instantiate the given resource.  Returns 0
	/// if the default ResourceOptions specified by the appropriate ResourceLoader
	/// should be used.
	virtual
	ResourceOptionsOP
	get_resource_options( ResourceDescription const & desc ) const = 0;

	/// @brief If the fallback configuration could not create a resource for the given
	/// resource description, then create an error message informing the user what options
	/// must be provided on the command line.
	virtual
	std::string
	could_not_create_resource_error_message( ResourceDescription const & desc ) const = 0;

};

} // namespace resource_manager
} // namespace basic

#endif // INCLUDED_basic_resource_manager_fallback_configuration_HH
