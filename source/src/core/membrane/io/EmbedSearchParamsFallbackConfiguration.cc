// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedSearchParamsFallbackConfiguration.cc
///
/// @brief      Embedding Search Parameters Fallback Configuration class  - Contains options for membrane search and score
/// @details    Membrane proteins in rosetta use the membrane scoring function and an MCM embedidng search
///             which can be tuned and adjusted using the following options.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedSearchParamsFallbackConfiguration_cc
#define INCLUDED_core_membrane_io_EmbedSearchParamsFallbackConfiguration_cc

// Unit Headers
#include <core/membrane/io/EmbedSearchParamsFallbackConfiguration.hh>
#include <core/membrane/io/EmbedSearchParamsFallbackConfigurationCreator.hh>

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

/// @brief Constructor
EmbedSearchParamsFallbackConfiguration::EmbedSearchParamsFallbackConfiguration() {}

/// @brief Allow no fallback
bool
EmbedSearchParamsFallbackConfiguration::fallback_specified( ResourceDescription const & ) const
{
    return false;
}

/// @brief Specify corresponding loader type
basic::resource_manager::LoaderType
EmbedSearchParamsFallbackConfiguration::get_resource_loader( ResourceDescription const & ) const
{
    return "EmbedSearchParams";
}

/// @brief Specify locator ID
basic::resource_manager::LocatorID
EmbedSearchParamsFallbackConfiguration::get_locator_id( ResourceDescription const & ) const
{
    return "";
}

/// @brief Specify corresponding resource options
basic::resource_manager::ResourceOptionsOP
EmbedSearchParamsFallbackConfiguration::get_resource_options( ResourceDescription const & ) const
{
    // use default options
    return NULL;
}

/// @brief Throw error message if someone tries to bypass .xml and go via commandline
std::string
EmbedSearchParamsFallbackConfiguration::could_not_create_resource_error_message( ResourceDescription const & ) const
{
    return "No fallback for embedding search parameters!";
}

/// @brief Creator Class - Create fallback configuration class
basic::resource_manager::FallbackConfigurationOP
EmbedSearchParamsFallbackConfigurationCreator::create_fallback_configuration() const
{
    return new EmbedSearchParamsFallbackConfiguration();
}

/// @brief Supply resource description to registrator
std::string
EmbedSearchParamsFallbackConfigurationCreator::resource_description() const
{
    return "EmbedSearchParams";
}

} // io
} // memrbane
} // core

#endif // INCLUDED_core_membrane_io_EmbedSearchParamsFallbackConfiguration_cc


