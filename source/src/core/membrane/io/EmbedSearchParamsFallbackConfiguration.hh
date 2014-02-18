// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedSearchParamsFallbackConfiguration.hh
///
/// @brief      Embedding Search Parameters Fallback Configuration class  - Contains options for membrane search and score
/// @details    Membrane proteins in rosetta use the membrane scoring function and an MCM embedidng search
///             which can be tuned and adjusted using the following options.
///
/// @note       Resource Manager Component
/// @note       This class is a dummy placeholder - not allowing a fallback for embedding search parameters*****
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_memrbane_io_EmbedSearchParamsFallbackConfiguration_hh
#define INCLUDED_core_memrbane_io_EmbedSearchParamsFallbackConfiguration_hh

// Unit Headers
#include <core/membrane/io/EmbedSearchParamsFallbackConfiguration.fwd.hh>
#include <basic/resource_manager/FallbackConfiguration.hh>
#include <basic/resource_manager/types.hh>

namespace core {
namespace membrane {
namespace io {
    
    /// @brief Embedding Search Parameters Fallback Configuration
    /// @details fallback if no resource for embedding search parameters is specified
    class EmbedSearchParamsFallbackConfiguration : public basic::resource_manager::FallbackConfiguration
    {
        
    public: // typedefs
        
        typedef basic::resource_manager::ResourceDescription ResourceDescription;
        
    public: // functions
        
        /// @brief Constructor
        EmbedSearchParamsFallbackConfiguration();
        
        /// @brief Specify fallback option on the commandline
        virtual
        bool
        fallback_specified( ResourceDescription const & desc ) const;
        
        /// @brief Specify corresponding resource loader for embedding search parameters
        virtual
        basic::resource_manager::LoaderType
        get_resource_loader( ResourceDescription const & desc ) const;
        
        /// @brief Specify locator id for resource (null)
        virtual
        basic::resource_manager::LocatorID
        get_locator_id( ResourceDescription const & desc ) const;
        
        /// @brief Specify corresponding resource options class for embed search parameters
        virtual
        basic::resource_manager::ResourceOptionsOP
        get_resource_options( ResourceDescription const & desc ) const;
        
        /// @brief Throw error message - fallback not allowed for embedding search parameters
        virtual
        std::string
        could_not_create_resource_error_message( ResourceDescription const & desc ) const;
        
    private:
        
        /// @brief Grab locator id (null)
        basic::resource_manager::LocatorID get_span_filename_from_options() const;
        
    }; // class EmbedSearchParamsFallbackConfiguraiton
    
} // io
} // membrane
} // core

#endif // INCLUDED_core_memrbane_io_EmbedSearchParamsFallbackConfiguration_hh

