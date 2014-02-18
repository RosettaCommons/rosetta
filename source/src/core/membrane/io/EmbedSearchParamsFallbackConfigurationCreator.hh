// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedSearchParamsFallbackConfigurationCreator.hh
///
/// @brief      Embedding Search Parameters Fallback Configuration class  - Contains options for membrane search and score
/// @details    Membrane proteins in rosetta use the membrane scoring function and an MCM embedidng search
///             which can be tuned and adjusted using the following options.
///
/// @note       Resource Manager Component
/// @note       This class is a dummy placeholder - not allowing a fallback for embedding search parameters*****
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedSearchParamsFallbackConfigurationCreator_hh
#define INCLUDED_core_membrane_io_EmbedSearchParamsFallbackConfigurationCreator_hh

//unit headers
#include <basic/resource_manager/FallbackConfigurationCreator.hh>

// package headers
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>

//C++ headers
#include <string>

namespace core {
namespace membrane {
namespace io {
    
    
    /// @brief Embedding Search Parameters - Fallback Configuration Registrator
    /// @details Register fallback config with embedding search params with core init
    class EmbedSearchParamsFallbackConfigurationCreator : public basic::resource_manager::FallbackConfigurationCreator
    {
    public:
        
        /// @brief Return fallback configuraiton class
        virtual
        basic::resource_manager::FallbackConfigurationOP
        create_fallback_configuration() const;
        
        /// @brief Return fallback configuration corresponding resource type
        virtual
        std::string resource_description() const;
        
    };
    
} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedSearchParamsFallbackConfigurationCreator_hh

