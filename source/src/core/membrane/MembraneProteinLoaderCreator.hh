// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/membrane/MembraneProteinLoaderCreator.hh
/// @brief  loader creator for membrane proteins
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_MembraneProteinLoaderCreator_hh
#define INCLUDED_core_membrane_MembraneProteinLoaderCreator_hh

// unit headers
#include <basic/resource_manager/ResourceLoaderCreator.hh>

namespace core {
namespace membrane {

    /// @brief Membrane Protein Loader Creator Class (RM)
    class MembraneProteinLoaderCreator : public basic::resource_manager::ResourceLoaderCreator
    {
    public:
        
        /// @brief Create Resource Loader for Membrane Proteins
        virtual
        basic::resource_manager::ResourceLoaderOP
        create_resource_loader() const;
        
        /// @brief Return loader type (MembraneProtein)
        virtual
        std::string loader_type() const;
        
    };
    
} // namespace membrane
} // namespace core

#endif // INCLUDED_core_membrane_MembraneProteinLoaderCreator_hh

