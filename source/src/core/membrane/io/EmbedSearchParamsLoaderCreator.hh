// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedSearchParamsLoaderCreator.hh
///
/// @brief      Embedding Search Parameters loader class  - Contains options for membrane search and score
/// @details    Membrane proteins in rosetta use the membrane scoring function and an MCM embedidng search
///             which can be tuned and adjusted using the following options.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedSearchParamsLoaderCreator_hh
#define INCLUDED_core_membrane_io_EmbedSearchParamsLoaderCreator_hh

// Unit Headers
#include <basic/resource_manager/ResourceLoaderCreator.hh>

namespace core {
namespace membrane {
namespace io {
            
    /// @brief Embed Search Parameters Loader Creator
    /// @details Register embedding search parameters loader with registrator
    class EmbedSearchParamsLoaderCreator : public basic::resource_manager::ResourceLoaderCreator {
        
    public:
        
        /// @brief Return resource loader (EmbedSearchParams)
        virtual
        basic::resource_manager::ResourceLoaderOP
        create_resource_loader() const;
        
        /// @brief Return resource loader type (EmbedSearchParams)
        virtual
        std::string
        loader_type() const;
        
    }; // class EmbedSearchParamsLoader
    
} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedSearchParamsLoaderCreator_hh

