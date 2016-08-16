// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   basic/resource_manager/fallback_configuration/NativeFallbackConfiguration.cc
/// @author Matthew O'Meara mattjomeara@gmail.com

// Unit Headers
#include <basic/resource_manager/fallback_configuration/NativeFallbackConfiguration.hh>
#include <basic/resource_manager/fallback_configuration/NativeFallbackConfigurationCreator.hh>

// basic headers
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

namespace basic {
namespace resource_manager {
namespace fallback_configuration {

using basic::resource_manager::LoaderType;
using basic::resource_manager::LocatorID;
using basic::resource_manager::LocatorTag;
using basic::resource_manager::ResourceDescription;
using basic::resource_manager::ResourceTag;
using basic::resource_manager::ResourceOptionsTag;

NativeFallbackConfiguration::NativeFallbackConfiguration()
{}

bool
NativeFallbackConfiguration::fallback_specified( ResourceDescription const & ) const
{
	using namespace basic::options;
	return option[ OptionKeys::in::file::native ].user();
}

basic::resource_manager::LoaderType
NativeFallbackConfiguration::get_resource_loader( ResourceDescription const & ) const
{
	return "PoseFromPDB";
}

basic::resource_manager::LocatorID
NativeFallbackConfiguration::get_locator_id( ResourceDescription const & ) const
{
	return get_native_filename_from_options();
}

basic::resource_manager::ResourceOptionsOP
NativeFallbackConfiguration::get_resource_options( ResourceDescription const & ) const
{
	// use the default options.
	return 0;
}

std::string
NativeFallbackConfiguration::could_not_create_resource_error_message( ResourceDescription const & ) const
{
	return "The NativeFallbackConfiguration requires that the flag '-in:file:native' be set on the command line.";
}

basic::resource_manager::LocatorID
NativeFallbackConfiguration::get_native_filename_from_options() const
{
	using namespace basic::options;
	return option[ OptionKeys::in::file::native ].value();
}

basic::resource_manager::FallbackConfigurationOP
NativeFallbackConfigurationCreator::create_fallback_configuration() const
{
	return basic::resource_manager::FallbackConfigurationOP( new NativeFallbackConfiguration );
}

std::string
NativeFallbackConfigurationCreator::resource_description() const
{
	return "native";
}

} // namespace
} // namespace
} // namespace
