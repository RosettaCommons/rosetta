// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedSearchParamsOptionsCreator.hh
///
/// @brief      Embedding Search Parameters Creator class  - Contains options for membrane search and score
/// @details    Membrane proteins in rosetta use the membrane scoring function and an MCM embedidng search
///             which can be tuned and adjusted using the following options.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedSearchParamsOptionsCreator_hh
#define INCLUDED_core_membrane_io_EmbedSearchParamsOptionsCreator_hh

//package headers
#include <basic/resource_manager/ResourceOptionsCreator.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

namespace core {
namespace membrane {
namespace io {
    
    /// @brief Embedding Search Parameters Options Creator
    /// @details Register Embed Search params options with core init registrator
    class EmbedSearchParamsOptionsCreator : public basic::resource_manager::ResourceOptionsCreator
    {
        
    public:
        
        /// @brief Return embed search params otpions type
        virtual std::string options_type() const;
        
        /// @brief Create Embedding search params options class
        virtual basic::resource_manager::ResourceOptionsOP create_options() const;
        
    }; // class EmbedSearchParamsOptionsCreator
    
} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedSearchParamsOptionsCreator_hh

