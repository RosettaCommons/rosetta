// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/looops/LoopsFileFallbackConfiguration.hh
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_protocols_loops_loops_file_fallback_configuration_HH
#define INCLUDED_protocols_loops_loops_file_fallback_configuration_HH

// Unit Header
#include <protocols/loops/LoopsFileFallbackConfiguration.fwd.hh>
#include <basic/resource_manager/FallbackConfiguration.hh>
#include <basic/resource_manager/types.hh>

namespace protocols {
namespace loops {

/// @brief %LoopsFileFallbackConfiguration provides instructions to the ResourceManager for loading loops_files in the
/// absence of a resource definition file.
/// @details The class will confirm that a fallback can be created and will provide the necessary information to the
/// ResourceManager to fully construct the %resource.
class LoopsFileFallbackConfiguration : public basic::resource_manager::FallbackConfiguration {
public:
	typedef basic::resource_manager::ResourceDescription ResourceDescription;
public:
	/// @brief Construct the %LoopsFileFallbackConfiguration.
	LoopsFileFallbackConfiguration();

	/// @brief Determine if the fallback configuration has been specified and return true or false.
	virtual
	bool
	fallback_specified( ResourceDescription const & desc ) const;

	/// @brief Return the type of loader that is required for this %resource.
	virtual
	basic::resource_manager::LoaderType
	get_resource_loader( ResourceDescription const & desc ) const;

	/// @brief Return the %locator_id that will be used to construct this %resource.
	virtual
	basic::resource_manager::LocatorID
	get_locator_id( ResourceDescription const & desc ) const;

	/// @brief Return an owning pointer to a ResourceOptions instance to configure this %resource.
	virtual
	basic::resource_manager::ResourceOptionsOP
	get_resource_options( ResourceDescription const & desc ) const;

	/// @brief Return a string that should be displayed if the %resource could not be created.
	virtual
	std::string
	could_not_create_resource_error_message( ResourceDescription const & desc ) const;

private:
	/// @brief Find and return the %locator_id for this %resource from the options system.
	basic::resource_manager::LocatorID get_loops_filename_from_options() const;

};

} // namespace loops
} // namespace protocols

#endif // INCLUDED_protocols_loops_loops_file_fallback_configuration_HH
