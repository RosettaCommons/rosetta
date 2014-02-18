// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/membrane/MembraneProteinFallbackConfigurationCreator.hh
/// @brief  Membrane protein fallback configuration creator
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneProteinFallbackConfiguration_creator_hh
#define INCLUDED_core_membrane_MembraneProteinFallbackConfiguration_creator_hh

//unit headers
#include <basic/resource_manager/FallbackConfigurationCreator.hh>

// package headers
#include <basic/resource_manager/FallbackConfiguration.fwd.hh>

//C++ headers
#include <string>

namespace core {
namespace membrane {

    /// @brief Membrane Protein Fallback Configuration Creator
    class MembraneProteinFallbackConfigurationCreator : public basic::resource_manager::FallbackConfigurationCreator
    {
    public:
        
        /// @brief Create a membrane protein fallback configuraiton - register with init.hh
        virtual
        basic::resource_manager::FallbackConfigurationOP
        create_fallback_configuration() const;
        
        /// @brief Return corresponding resource description
        virtual
        std::string resource_description() const;
        
    };
        
} // namespace membrane
} // namespace core

#endif // INCLUDED_core_membrane_MembraneProteinFallbackConfiguration_creator_hh

