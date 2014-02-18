// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedSearchParamsLoader.cc
///
/// @brief      Embedding Search Parameters loader class  - Contains options for membrane search and score
/// @details    Membrane proteins in rosetta use the membrane scoring function and an MCM embedidng search
///             which can be tuned and adjusted using the following options.
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedSearchParamsLoader_cc
#define INCLUDED_core_membrane_io_EmbedSearchParamsLoader_cc

// Unit headers
#include <core/membrane/io/EmbedSearchParamsLoader.hh>
#include <core/membrane/io/EmbedSearchParamsLoaderCreator.hh>

// Project Headers
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/definitions_util.hh>

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
    
    /// @brief Constructor
    EmbedSearchParamsLoader::EmbedSearchParamsLoader() {}

    /// @brief Destructor
    EmbedSearchParamsLoader::~EmbedSearchParamsLoader() {}
    
    /// @brief Create embedding search data object from resource
    utility::pointer::ReferenceCountOP
    EmbedSearchParamsLoader::create_resource(
                                    basic::resource_manager::ResourceOptions const & options,
                                    basic::resource_manager::LocatorID const &,
                                    std::istream &
                                    ) const
    {
        using namespace core::membrane::util;
        using namespace core::membrane::io;
        
        // Embedding Search Options
        if ( ! dynamic_cast< EmbedSearchParamsOptions const * > ( &options ) ) {
            throw utility::excn::EXCN_Msg_Exception("SpanFileLoader excpected to be given an EmbedSearchParamsOptions but was given a non-EmbedSearchParamsOptions object of type " + options.type() + "' which has the name '" + options.name() + "'.");
        }
        EmbedSearchParamsOptions const & opts = static_cast< EmbedSearchParamsOptions const & > ( options );

        
        // Load an embed search object from file
        EmbedSearchParamsIO esfio;
        EmbedSearchParamsOP params = esfio.get_embed_params_from_file( opts );
        
        // Return search params data
        return params;
    }
    
    /// @brief Import a default options object
    basic::resource_manager::ResourceOptionsOP
    EmbedSearchParamsLoader::default_options() const
    {
        using namespace core::membrane::io;
        
        return new EmbedSearchParamsOptions();
    }
    
    /// @brief Returns a new resource loader
    basic::resource_manager::ResourceLoaderOP
    EmbedSearchParamsLoaderCreator::create_resource_loader() const
    {
        using namespace core::membrane::io;
        
        return new EmbedSearchParamsLoader();
    }
    
    /// @brief Embed Definition loader type
    std::string
    EmbedSearchParamsLoaderCreator::loader_type() const
    {
        using namespace core::membrane::io;
        
        return "EmbedSearchParams";
    }
    
} // io
} // membrane
} // core 

#endif // INCLUDED_core_membrane_io_EmbedSearchParamsLoader_cc

