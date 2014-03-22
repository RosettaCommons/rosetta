// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedSearchParamsLoader.hh
///
/// @brief      Embedding Search Parameters loader class  - Contains options for membrane search and score
/// @details    Membrane proteins in rosetta use the membrane scoring function and an MCM embedidng search
///             which can be tuned and adjusted using the following options.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedSearchParamsLoader_hh
#define INCLUDED_core_membrane_io_EmbedSearchParamsLoader_hh

// Unit headers
#include <core/membrane/io/EmbedSearchParamsLoader.fwd.hh>

// Project Headers
#include <basic/resource_manager/ResourceLoader.hh>
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/resource_manager/types.hh>

#include <core/conformation/membrane/definitions.hh>
#include <core/conformation/membrane/definitions_util.hh>

#include <core/membrane/io/EmbedSearchParamsOptions.hh>
#include <core/membrane/io/EmbedSearchParamsIO.hh>

// Utility Headers
#include <utility/excn/Exceptions.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <cstdlib>

namespace core {
namespace membrane {
namespace io {

/// @brief Reosurce Loader for Embedding Search Parameters
/// @details Loads embedding search and score parameters as a resource
class EmbedSearchParamsLoader : public basic::resource_manager::ResourceLoader {
        
public:
        
    /// @brief Constructor
    EmbedSearchParamsLoader();
    
    /// @brief Destructor
    virtual ~EmbedSearchParamsLoader();
    
    /// @brief Returns an OP to embedding search options
    virtual
    utility::pointer::ReferenceCountOP
    create_resource(
                    basic::resource_manager::ResourceOptions const & options,
                    basic::resource_manager::LocatorID const & locator_id,
                    std::istream & istream
                    ) const;
    
    /// @brief Default options
    virtual
    basic::resource_manager::ResourceOptionsOP
    default_options() const;
    
}; // class EmbedSearchParamsLoader
    
} // io
} // membrane
} // core 

#endif // INCLUDED_core_membrane_io_EmbedSearchParamsLoader_hh

