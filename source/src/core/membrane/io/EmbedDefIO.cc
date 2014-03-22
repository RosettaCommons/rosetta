// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedDefIO.cc
///
/// @brief      Read in data from files for membrane protein chain embeddings
/// @details    Protein embedding is defined by a normal and center vector positioned
///             with respect to the membrane as well as a specified depth
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedDefIO_cc
#define INCLUDED_core_membrane_io_EmbedDefIO_cc

// Unit Headers
#include <core/membrane/io/EmbedDefIO.hh>

// Project Headers
#include <core/conformation/membrane/definitions.hh>
#include <core/conformation/membrane/definitions_util.hh>
#include <core/conformation/membrane/Exceptions.hh>

// Package Headers
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <algorithm>
#include <string>
#include <cstdlib>
#include <cmath>

using basic::Error;
using basic::Warning;

using std::stringstream;
using std::string;

static basic::Tracer TR( "core.membrane.io.EmbedDefIO" );

using namespace core::membrane;
using namespace core::conformation::membrane;

namespace core {
namespace membrane {
namespace io {

    /// @brief Constructor
    EmbedDefIO::EmbedDefIO() :
        utility::pointer::ReferenceCount()
    {}

    /// @brief Destructor
    EmbedDefIO::~EmbedDefIO()
    {}

    /// @brief Read and store data from embedding definition file
    /// @param embedding definition object, embedfile
    void EmbedDefIO::read_embed_data(
                                    EmbedConfigInfoOP embedding,
                                    std::string embedfile
                                    )
    {

        using namespace core;

        TR << "Initializing embedding definition info using " << embedfile << std::endl;

        // Create an izstream for the file
        std::string line;
        utility::io::izstream stream (embedfile);

        // Open file
        stream.open(embedfile);

        if (stream) {

            // Initialize some Vars
            std::string param;
            Real x; Real y; Real z;
            std::string tag;

            // Get Pose Description Line
            getline(stream, line);

            // Read in String Stream Normal Line
            getline(stream, line);
            std::istringstream l(line);
            l >> param;
            l >> x; l >> y; l >> z;
            l >> tag;

            // Set the values in the embedding configuration
            embedding->normal.assign(x, y, z);
            embedding->normal_tag = tag;

            // Read in String Stream center Line
            getline(stream, line);
            std::istringstream m(line);
            m >> param;
            m >> x; m >> y; m >> z;
            m >> tag;

            // Set the values in the embedding configuration
            embedding->center.assign(x, y, z);
            embedding->center_tag = tag;

            // Read in the String Stream for Depth
            getline(stream, line);
            std::istringstream n(line);
            n >> param;
            n >> x;

            // Set value in embedding configuration
            embedding->depth = x;

        } else {
            throw new EXCN_Illegal_Arguments("Embedding data file not found!");
        }

    }

    /// @brief Get Embedding Definition from user specified file
    EmbedConfigInfoOP
    EmbedDefIO::get_embedding_from_file( std::string embedfile )
    {

        using namespace core::membrane::io;

        if ( embedfile.compare("") == 0 ) {
            throw new EXCN_Illegal_Arguments("Illegal embedding definition file locator");
        }

        /// Create new object
        EmbedConfigInfoOP embedding = init_embedConfigInfo();

        // read in Lips data
        read_embed_data( embedding, embedfile );

        // Done!
        return embedding;
    }

} // geometry
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedDefIO_cc

