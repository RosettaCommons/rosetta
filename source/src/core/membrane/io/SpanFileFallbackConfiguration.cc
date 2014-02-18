// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       SpanFileFallbackConfiguration.cc
///
/// @brief      Fallback configuration for span file loader/options - generating per-chain membrane spanning data
/// @details    Generates membrane spanning topology data from Octopus topology
///             prediction information. Topology can be generated at http://octopus.cbr.su.se/
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_SpanFileFallbackConfiguration_cc
#define INCLUDED_core_membrane_io_SpanFileFallbackConfiguration_cc

// Unit Headers
#include <core/membrane/io/SpanFileFallbackConfiguration.hh>
#include <core/membrane/io/SpanFileFallbackConfigurationCreator.hh>

// Platform Headers
#include <core/types.hh>
#include <utility/vector1.hh>

// Basic Headers
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Nameric Headers
#include <numeric/random/random.hh>

// Utility Headers
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
    
/// @brief Fallback Configuration Constructor
SpanFileFallbackConfiguration::SpanFileFallbackConfiguration() {}

/// @brief Return fallback option -in::file::spanfile
bool
SpanFileFallbackConfiguration::fallback_specified( ResourceDescription const & ) const
{
	return basic::options::option[ basic::options::OptionKeys::in::file::spanfile ].user();
}

/// @brief Return corresponding loader type - SpanFile
basic::resource_manager::LoaderType
SpanFileFallbackConfiguration::get_resource_loader( ResourceDescription const & ) const
{
	return "SpanFile";
}

/// @brief Return locator id of source
basic::resource_manager::LocatorID
SpanFileFallbackConfiguration::get_locator_id( ResourceDescription const & ) const
{
	return get_span_filename_from_options();
}

/// @brief Get corresponding options class for Span Files
basic::resource_manager::ResourceOptionsOP
SpanFileFallbackConfiguration::get_resource_options( ResourceDescription const & ) const
{
	// use default options
	return NULL;
}

/// @brief Throw error if fallback or options not specified
std::string
SpanFileFallbackConfiguration::could_not_create_resource_error_message( ResourceDescription const & ) const
{
	return "The SpanFileFallbackConfiguration requires the flag -in::file::spanfile to be set on the command line";
}

/// @brief Get span file from options
basic::resource_manager::LocatorID
SpanFileFallbackConfiguration::get_span_filename_from_options() const
{

	std::string spanfile = basic::options::option[ basic::options::OptionKeys::in::file::spanfile ]();
	if ( spanfile.compare("") == 0 )
	{
		throw utility::excn::EXCN_Msg_Exception("The fallback SpanFile resource option has no span file associated with it! Was the option omitted from the command line?");
	}

	return spanfile;
}

/// @brief Creator Class - Return Fallback Configuraiton Class
basic::resource_manager::FallbackConfigurationOP
SpanFileFallbackConfigurationCreator::create_fallback_configuration() const
{
	return new SpanFileFallbackConfiguration();
}
    
/// @brief Creator Class - Return corresponding loader/resource type
std::string
SpanFileFallbackConfigurationCreator::resource_description() const
{
	return "SpanFile";
}

} // io
} // memrbane
} // core

#endif // INCLUDED_core_membrane_io_SpanFileFallbackConfiguration_cc
