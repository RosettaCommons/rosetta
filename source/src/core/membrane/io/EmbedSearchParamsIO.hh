// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedSearchParamsIO.cc
///
/// @brief      Embedding Search Parameters IO class  - Contains options for membrane search and score
/// @details    Membrane proteins in rosetta use the membrane scoring function and an MCM embedidng search
///             which can be tuned and adjusted using the following options.
///
/// @note       Resource Manager Component
/// @note       Container class - reads straight from options class
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedSearchParamsIO_hh
#define INCLUDED_core_membrane_io_EmbedSearchParamsIO_hh

// Unit Headers
#include <core/membrane/io/EmbedSearchParamsIO.fwd.hh>

// Project Headers
#include <core/conformation/membrane/definitions.hh>
#include <core/conformation/membrane/definitions_util.hh>
#include <core/conformation/membrane/Exceptions.hh>

#include <core/membrane/io/EmbedSearchParamsOptions.hh>

// Package Headers
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <utility/tag/Tag.hh>

// C++ Headers
#include <algorithm>
#include <string>
#include <cstdlib>
#include <cmath>

using namespace core::conformation::membrane;

namespace core {
namespace membrane {
namespace io {
            
    /// @brief Embeding Search Params IO Class
    /// @details Read data straight from options class and store as a resource
    class EmbedSearchParamsIO : public utility::pointer::ReferenceCount
    {
        
    public:
        
        /// @brief Constructor
        EmbedSearchParamsIO();
        
        /// @brief Destructor
        ~EmbedSearchParamsIO();
        
        
        /// @brief Main IO Function - reads data straight from options class
        EmbedSearchParamsOP
        get_embed_params_from_file( core::membrane::io::EmbedSearchParamsOptions const & opts );
        
    };

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedSearchParamsIO_hh

