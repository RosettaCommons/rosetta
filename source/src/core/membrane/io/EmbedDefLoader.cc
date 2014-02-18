// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedDefLoader.cc
///
/// @brief      Loader class for membrane protein chain embeddings
/// @details    Protein embedding is defined by a normal and center vector positioned
///             with respect to the membrane as well as a specified depth
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedDefLoader_cc
#define INCLUDED_core_membrane_io_EmbedDefLoader_cc

// Unit headers
#include <core/membrane/io/EmbedDefLoader.hh>
#include <core/membrane/io/EmbedDefLoaderCreator.hh>

// Project Headers
#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/definitions_util.hh>

#include <core/membrane/io/EmbedDefOptions.hh>
#include <core/membrane/io/EmbedDefIO.hh>

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
    EmbedDefLoader::EmbedDefLoader() {}
    
    /// @brief Destructor
    EmbedDefLoader::~EmbedDefLoader() {}

    /// @brief Create embedding search data object from resource
    utility::pointer::ReferenceCountOP
    EmbedDefLoader::create_resource(
        basic::resource_manager::ResourceOptions const &,
        basic::resource_manager::LocatorID const & locator_id,
        std::istream &
        ) const
    {
        using namespace core::membrane::util;
        using namespace core::membrane::io;

        // Load an embed search object from file
        EmbedDefIO edfio;
        EmbedConfigInfoOP embedding = edfio.get_embedding_from_file( locator_id );

        // Return initialized embedding data object
        return embedding;
    }

    /// @brief Import a default embedding options object
    basic::resource_manager::ResourceOptionsOP
    EmbedDefLoader::default_options() const
    {
        using namespace core::membrane::io;
        
        return new EmbedDefOptions();
    }

    /// @brief Returns a new resource loader
    basic::resource_manager::ResourceLoaderOP
    EmbedDefLoaderCreator::create_resource_loader() const
    {
        using namespace core::membrane::io;
        
        return new EmbedDefLoader();
    }

    /// @brief Embed Definition loader type
    std::string
    EmbedDefLoaderCreator::loader_type() const
    {
        using namespace core::membrane::io;
        
        return "EmbedDef";
    }

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedDefLoader_cc

