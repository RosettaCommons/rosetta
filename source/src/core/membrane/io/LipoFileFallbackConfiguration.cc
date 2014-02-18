// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       LipoFileFallbackConfiguration.cc
///
/// @brief      Fallback Configuration class - load object that stores Lipid Exposure Data from Membrane Spanning Topology
/// @details    Lipid exposure data calculated using run_lips.pl from octopus spanning
///             topology data. Requires blast.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_LipoFileFallbackConfiguration_cc
#define INCLUDED_core_membrane_io_LipoFileFallbackConfiguration_cc

// Unit Headers
#include <core/membrane/io/LipoFileFallbackConfigurationCreator.hh>
#include <core/membrane/io/LipoFileFallbackConfiguration.hh>

// Platform headers
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/definitions_util.hh>

// Basic Headers
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Numeric Header
#include <numeric/random/random.hh>

// Utility Header
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>

// C++ Headers
#include <string>
#include <map>

namespace core {
namespace membrane {
namespace io {

using basic::resource_manager::LoaderType;
using basic::resource_manager::LocatorID;
using basic::resource_manager::LocatorTag;
using basic::resource_manager::ResourceDescription;
using basic::resource_manager::ResourceTag;
using basic::resource_manager::ResourceOptionsTag;

/// @brief Constructor
LipoFileFallbackConfiguration::LipoFileFallbackConfiguration()
{}

/// @brief Option to Fallback from (-in:file:lipofile)
bool
LipoFileFallbackConfiguration::fallback_specified( ResourceDescription const & ) const
{
	return basic::options::option[ basic::options::OptionKeys::in::file::lipofile ].user();
}

/// @brief Return loader type
basic::resource_manager::LoaderType
LipoFileFallbackConfiguration::get_resource_loader( ResourceDescription const & ) const
{
	return "LipoFile";
}

/// @brief Return locator id for lipo file
basic::resource_manager::LocatorID
LipoFileFallbackConfiguration::get_locator_id( ResourceDescription const & ) const
{
	return basic::options::option[ basic::options::OptionKeys::in::file::lipofile ].value();
}

/// @brief Return lipo file options
basic::resource_manager::ResourceOptionsOP
LipoFileFallbackConfiguration::get_resource_options( ResourceDescription const & ) const
{
	// use the default surface_vector options.
	return NULL;
}

/// @brief Throw error if no resource specified
std::string
LipoFileFallbackConfiguration::could_not_create_resource_error_message( ResourceDescription const & ) const
{
	return "The LipoFileFallbackConfiguration requires that the flag '-in:file:lipofile' be set on the command line.";
}

/// @brief Return lipo file fallback configuration class
basic::resource_manager::FallbackConfigurationOP
LipoFileFallbackConfigurationCreator::create_fallback_configuration() const
{
	return new LipoFileFallbackConfiguration;
}

/// @brief Return fallback configuration type
std::string
LipoFileFallbackConfigurationCreator::resource_description() const
{
	return "lipofile";
}

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_LipoFileFallbackConfiguration_cc

