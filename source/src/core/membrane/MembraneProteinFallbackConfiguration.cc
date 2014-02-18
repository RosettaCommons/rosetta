// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/membrane/MembraneProteinFallbackConfiguration.cc
/// @brief  Membrane protein fallback configuration
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneProteinFallbackConfiguration_cc
#define INCLUDED_core_membrane_MembraneProteinFallbackConfiguration_cc

// Unit Headers
#include <core/membrane/MembraneProteinFallbackConfiguration.hh>
#include <core/membrane/MembraneProteinFallbackConfigurationCreator.hh>

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

//C++ Headers
#include <cstdlib>
#include <string>
#include <map>

namespace core {
namespace membrane {
    
    using basic::resource_manager::LoaderType;
    using basic::resource_manager::LocatorID;
    using basic::resource_manager::LocatorTag;
    using basic::resource_manager::ResourceDescription;
    using basic::resource_manager::ResourceTag;
    using basic::resource_manager::ResourceOptionsTag;

    
    /// @brief Constructor
    MembraneProteinFallbackConfiguration::MembraneProteinFallbackConfiguration()
    {}
    
    /// @brief Specify a fallback option and result if applicable
    bool
    MembraneProteinFallbackConfiguration::fallback_specified( ResourceDescription const & ) const
    {
        return basic::options::option[ basic::options::OptionKeys::in::membrane ].user();
    }
    
    /// @brief Return applicable loader type for the resource
    basic::resource_manager::LoaderType
    MembraneProteinFallbackConfiguration::get_resource_loader( ResourceDescription const & ) const
    {
        return "MembraneProtein";
    }
    
    /// @brief Return locator ID from which to load the new resource
    basic::resource_manager::LocatorID
    MembraneProteinFallbackConfiguration::get_locator_id( ResourceDescription const & ) const
    {
        return "";
    }
    
    /// @brief Return resource options class for membrane proteins
    basic::resource_manager::ResourceOptionsOP
    MembraneProteinFallbackConfiguration::get_resource_options( ResourceDescription const & ) const
    {
        // use the default loops file options.
        return 0;
    }
    
    /// @brief Returns an error message if no fallback or loader specified
    std::string
    MembraneProteinFallbackConfiguration::could_not_create_resource_error_message( ResourceDescription const & ) const
    {
        std::string const msg = "The MembraneProteinFallbackConfiguration requires that the flag '-in:membrane' be set on the command line.";
    
        return msg;
    }
    
    /// @brief Get options from options class
    bool
    MembraneProteinFallbackConfiguration::get_membrane_opts_from_options() const
    {
        bool const use_membrane = basic::options::option[ basic::options::OptionKeys::in::membrane ]();
        return use_membrane;
    }
    
    /// @brief Create fallback configuraiton class for membrane protein
    basic::resource_manager::FallbackConfigurationOP
    MembraneProteinFallbackConfigurationCreator::create_fallback_configuration() const
    {
        return new MembraneProteinFallbackConfiguration;
    }
    
    /// @brief Return resource description applicable
    std::string
    MembraneProteinFallbackConfigurationCreator::resource_description() const
    {
        return "MembraneProtein";
    }

        
} // namespace membrane
} // namespace core

#endif // INCLUDED_core_membrane_MembraneProteinFallbackConfiguration_cc

