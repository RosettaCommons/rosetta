// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/FallbackConfiguration.hh
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_basic_resource_manager_fallback_configuration_HH
#define INCLUDED_basic_resource_manager_fallback_configuration_HH

// Unit Headers
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>

// Package headers
#include <basic/resource_manager/ResourceOptions.fwd.hh>
#include <basic/resource_manager/types.hh>


namespace basic {
namespace resource_manager {

/// @brief The %FallbackConfiguration class describes how a resource should
/// be created if the ResourceManager is not being used, i.e., if resources
/// have been specified through flags on the command line instead of through
/// a resource-definition file.
///
/// A protocol will still query the ResourceManager
/// for its desired resource in such a case (by its "resource_description")
/// and if no resource has been provided for that description, then the
/// ResourceManager will ask the FallbackConfigurationFactory for the
/// %FallbackConfiguration that has been registered for that resource
/// description (this registration logic is handled by the
/// FallbackConfigurationCreator) and ask it to provide the information
/// necessary to create the appropriate resource from the command line.
///
/// For example, a user can define the loop they are interested in modeling to
/// the loop-modeling protocol via the LoopsFileData resource.  The loop modeling
/// protocol once was but is no longer responsible for opening a text file
/// and creating a LoopsFileData object from that text file; instead it asks the
/// ResourceManager directly for the resource.  For the sake of backwards
/// compatibility, the ResourceManager must have some way of creating a
/// LoopsFileData object from the command line.  If someone is using the loop modeling
/// application, but does not provide a resource definition, then the
/// LoopsFileFallbackConfiguration (derived from %FallbackConfiguration)
/// will be called upon by the ResourceManager for three things:
/// 1) the name of the ResourceLoader to use to create the resource ("LoopsFile"),
/// 2) the locator_id for the Resource (which it goes to the command line to find), and
/// 3) a ResourceOptionsOP which is either null (meaning "use the default") or is in
/// some way initialized with options from the command line.  With these three bits, the
/// ResourceManager will go the appropriate ResourceLoader and request the construction
/// of the Resource as usual, opening the file using the default ResourceLocator
/// (for the JD2ResourceManager, this is the FileSystemResourceLocator).
class FallbackConfiguration : public utility::pointer::ReferenceCount {
public:
	/// @brief Default constructor for the FallbackConfiguration initializes nothing
	FallbackConfiguration() {}

	/// @brief Destructor for the FallbackConfiguration does nothing
	virtual ~FallbackConfiguration();

	/// @brief Return true if a fallback been provided for the given resource description.
	/// It is possible that the fallback configuration would be unable to deliver
	/// a requested resource, e.g. if the appropriate command line option has not been provided.
	virtual
	bool
	fallback_specified( ResourceDescription const & desc ) const = 0;

	/// @brief Return the name of the resource loader that should be used to instantiate the
	/// desired resource
	virtual
	LoaderType
	get_resource_loader( ResourceDescription const & desc ) const = 0;

	/// @brief Return the locator id that the resource locator should use to
	/// find the data used to construct the desired resource - e.g. the file name.
	virtual
	LocatorID
	get_locator_id( ResourceDescription const & desc ) const = 0;

	/// @brief Return a pointer to the ResourceOptions object that the
	/// ResourceLoader will use to instantiate the given resource.  Return 0
	/// if the default ResourceOptions specified by the appropriate ResourceLoader
	/// should be used.
	virtual
	ResourceOptionsOP
	get_resource_options( ResourceDescription const & desc ) const = 0;

	/// @brief If the %FallbackConfiguration has not been provided the appropriate set
	/// of command line flags needed to construct the Resource given the
	/// resource description, then return an error message informing the user
	/// what options must be provided on the command line.
	virtual
	std::string
	could_not_create_resource_error_message( ResourceDescription const & desc ) const = 0;

};

} // namespace resource_manager
} // namespace basic

#endif // INCLUDED_basic_resource_manager_fallback_configuration_HH
