// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/surface_docking/SurfaceVectorFallbackConfiguration.cc
/// @author Michael Pacella mpacella88@gmail.com

// Unit Headers
#include <protocols/surface_docking/SurfaceVectorFallbackConfiguration.hh>
#include <protocols/surface_docking/SurfaceVectorFallbackConfigurationCreator.hh>


// Platform Headers
#include <core/types.hh>
#include <utility/vector1.hh>

// basic headers
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// numeric headers
#include <numeric/random/random.hh>

//utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>

//C++ Headers
#include <string>
#include <map>


namespace protocols {
namespace surface_docking {

using basic::resource_manager::LoaderType;
using basic::resource_manager::LocatorID;
using basic::resource_manager::LocatorTag;
using basic::resource_manager::ResourceDescription;
using basic::resource_manager::ResourceTag;
using basic::resource_manager::ResourceOptionsTag;



SurfaceVectorFallbackConfiguration::SurfaceVectorFallbackConfiguration()
{}

bool
SurfaceVectorFallbackConfiguration::fallback_specified( ResourceDescription const & ) const
{
	return basic::options::option[ basic::options::OptionKeys::in::file::surface_vectors ].user();
}

basic::resource_manager::LoaderType
SurfaceVectorFallbackConfiguration::get_resource_loader( ResourceDescription const & ) const
{
	return "SurfaceVector";
}

basic::resource_manager::LocatorID
SurfaceVectorFallbackConfiguration::get_locator_id( ResourceDescription const & ) const
{
	return basic::options::option[ basic::options::OptionKeys::in::file::surface_vectors ].value();
}

basic::resource_manager::ResourceOptionsOP
SurfaceVectorFallbackConfiguration::get_resource_options( ResourceDescription const & ) const
{
	// use the default surface_vector options.
	return NULL;
}

std::string
SurfaceVectorFallbackConfiguration::could_not_create_resource_error_message( ResourceDescription const & ) const
{
	return "The SurfaceVectorFallbackConfiguration requires that the flag '-in:file:surface_vectors' be set on the command line.";
}


basic::resource_manager::FallbackConfigurationOP
SurfaceVectorFallbackConfigurationCreator::create_fallback_configuration() const
{
	return new SurfaceVectorFallbackConfiguration;
}

std::string
SurfaceVectorFallbackConfigurationCreator::resource_description() const
{
	return "surface_vectors";
}

} // namespace surface_docking
} // namespace protocols
