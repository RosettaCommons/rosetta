// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedDefFallbackConfiguration.cc
///
/// @brief      Fallback configuration for embed def - membrane protein chain embeddings
/// @details    Protein embedding is defined by a normal and center vector positioned
///             with respect to the membrane as well as a specified depth
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedDefFallbackConfiguration_cc
#define INCLUDED_core_membrane_io_EmbedDefFallbackConfiguration_cc

// Unit Headers
#include <core/membrane/io/EmbedDefFallbackConfiguration.hh>
#include <core/membrane/io/EmbedDefFallbackConfigurationCreator.hh>

// Platform Headers
#include <core/types.hh>
#include <core/conformation/membrane/definitions.hh>

// Basic Headers
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Numeric Headers
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

/// @brief Constructor
EmbedDefFallbackConfiguration::EmbedDefFallbackConfiguration() {}

/// @brief Destructor
EmbedDefFallbackConfiguration::~EmbedDefFallbackConfiguration() {}

/// @brief Specify Fallback for embedding definition resource
/// @details fallback on option group key (-in:file:embedfile)
bool
EmbedDefFallbackConfiguration::fallback_specified( ResourceDescription const & ) const
{

	return basic::options::option[ basic::options::OptionKeys::in::file::embedfile ].user();

}

/// @brief Get resource loader type
basic::resource_manager::LoaderType
EmbedDefFallbackConfiguration::get_resource_loader( ResourceDescription const & ) const
{
	return "EmbedDef";
}

/// @brief Return locator id of specified commandline file
basic::resource_manager::LocatorID
EmbedDefFallbackConfiguration::get_locator_id( ResourceDescription const & ) const
{
	return get_embedfile_from_options();
}

/// @brief Get corresponding resource options for embedding definitions
basic::resource_manager::ResourceOptionsOP
EmbedDefFallbackConfiguration::get_resource_options( ResourceDescription const & ) const
{
	// default options
	return 0;
}

/// @brief Throw error if no resource specified
std::string
EmbedDefFallbackConfiguration::could_not_create_resource_error_message( ResourceDescription const & ) const
{
	return "The Embedding definition fallback configuration requires the flag -membrane:embed_def to be set on the command line";
}

/// @brief Grab Embedding Definition file locator id
basic::resource_manager::LocatorID
EmbedDefFallbackConfiguration::get_embedfile_from_options() const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// Not doing assignment here - at least for now
	if ( !option[ in::file::embedfile ].user() ) {
		throw utility::excn::EXCN_Msg_Exception("The fallback embedding definition requires you to actually specify values!" );
	}

	// Get the file
	std::string embedfile = option[ in::file::embedfile ]();

	// Return the locator id to the embed def file
	return embedfile;
}

/// @brief Creator class - Return new embed def fallback configuration
basic::resource_manager::FallbackConfigurationOP
EmbedDefFallbackConfigurationCreator::create_fallback_configuration() const
{
	return new EmbedDefFallbackConfiguration;
}

/// @brief Creator class - return fallback configuration type to registrator
std::string
EmbedDefFallbackConfigurationCreator::resource_description() const
{
	return "EmbedDef";
}

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedDefFallbackConfiguration_cc

