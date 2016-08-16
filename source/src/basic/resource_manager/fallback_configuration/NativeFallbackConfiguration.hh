// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/fallback_configurations/NativeFallbackConfiguration.hh
/// @author Matthew O'Meara mattjomeara@gmail.com

#ifndef INCLUDED_basic_resource_manager_fallback_configurations_NativeFallbackConfiguration_HH
#define INCLUDED_basic_resource_manager_fallback_configurations_NativeFallbackConfiguration_HH

// Unit Header
#include <basic/resource_manager/fallback_configuration/NativeFallbackConfiguration.fwd.hh>
#include <basic/resource_manager/FallbackConfiguration.hh>
#include <basic/resource_manager/types.hh>

namespace basic {
namespace resource_manager {
namespace fallback_configuration {

class NativeFallbackConfiguration : public basic::resource_manager::FallbackConfiguration {
public:
	typedef basic::resource_manager::ResourceDescription ResourceDescription;
public:
	NativeFallbackConfiguration();

	virtual
	bool
	fallback_specified( ResourceDescription const & desc ) const;

	virtual
	basic::resource_manager::LoaderType
	get_resource_loader( ResourceDescription const & desc ) const;

	virtual
	basic::resource_manager::LocatorID
	get_locator_id( ResourceDescription const & desc ) const;

	virtual
	basic::resource_manager::ResourceOptionsOP
	get_resource_options( ResourceDescription const & desc ) const;

	virtual
	std::string
	could_not_create_resource_error_message( ResourceDescription const & desc ) const;

private:

	basic::resource_manager::LocatorID get_native_filename_from_options() const;

};

} // namespace
} // namespace
} // namespace

#endif // include guard
