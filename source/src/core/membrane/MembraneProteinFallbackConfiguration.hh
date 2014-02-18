// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/membrane/MembraneProteinFallbackConfiguration.hh
/// @brief  Membrane protein fallback configuration
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneProteinFallbackConfiguration_hh
#define INCLUDED_core_membrane_MembraneProteinFallbackConfiguration_hh

// Unit Header
#include <core/membrane/MembraneProteinFallbackConfiguration.fwd.hh>
#include <basic/resource_manager/FallbackConfiguration.hh>
#include <basic/resource_manager/types.hh>

namespace core {
namespace membrane {
    
    /// @brief Class for membrane protein fallback configuration
    class MembraneProteinFallbackConfiguration : public basic::resource_manager::FallbackConfiguration {
    public:
        
        typedef basic::resource_manager::ResourceDescription ResourceDescription;
    
    public:
        
        /// @brief Constructor
        MembraneProteinFallbackConfiguration();
        
        /// @brief Specify if a fallback was given for this resource
        virtual
        bool
        fallback_specified( ResourceDescription const & desc ) const;
        
        /// @brief Specify the corresponding loader type
        virtual
        basic::resource_manager::LoaderType
        get_resource_loader( ResourceDescription const & desc ) const;
        
        /// @brief Return appropriate locator id for teh resource given in the fallback option
        virtual
        basic::resource_manager::LocatorID
        get_locator_id( ResourceDescription const & desc ) const;
        
        /// @brief Provide a resource description
        virtual
        basic::resource_manager::ResourceOptionsOP
        get_resource_options( ResourceDescription const & desc ) const;
        
        /// @brief Return error message if no fallback or resource description specified
        virtual
        std::string
        could_not_create_resource_error_message( ResourceDescription const & desc ) const;
        
    private:
        
        /// @brief Get resource from options
        bool get_membrane_opts_from_options() const;
        
    };
    

} // namespace membrane
} // namespace core

#endif // INCLUDED_core_membrane_MembraneProteinFallbackConfiguration_hh

