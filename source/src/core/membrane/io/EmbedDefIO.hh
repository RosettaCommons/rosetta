// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       EmbedDefIO.hh
///
/// @brief      Read in data from files for membrane protein chain embeddings
/// @details    Protein embedding is defined by a normal and center vector positioned
///             with respect to the membrane as well as a specified depth
///
/// @note       Resource Manager Component
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_io_EmbedDefIO_hh
#define INCLUDED_core_membrane_io_EmbedDefIO_hh

// Unit Headers
#include <core/membrane/io/EmbedDefIO.fwd.hh>

// Project Headers
#include <core/conformation/membrane/definitions.hh>
#include <core/conformation/membrane/definitions_util.hh>
#include <core/conformation/membrane/Exceptions.hh>

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

/// @brief Embedding Definition File Reader
/// @details Reads a user-specified .emebed file containing starting embedding coordinates
class EmbedDefIO : public utility::pointer::ReferenceCount
{

public:

    /// @brief Constructor
    EmbedDefIO();
    
    /// @brief Destructor
    ~EmbedDefIO();
    
 
    /// @brief Read Embedding Definition from specified file
    EmbedConfigInfoOP
    get_embedding_from_file( std::string embedfile );
    
private:
    
    /// @brief Read and store data from lips file
    /// @param Lipid Info object, lipid info file, and pose
    void
    read_embed_data(
                                     EmbedConfigInfoOP embedding,
                                     std::string embedfile
                                     );

    
};

} // io
} // membrane
} // core

#endif // INCLUDED_core_membrane_io_EmbedDefIO_hh

