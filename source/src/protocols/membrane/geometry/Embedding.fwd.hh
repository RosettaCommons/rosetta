// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/membrane/geometry/Embedding.hh
///
/// @brief      Methods for Computing Membrane Embeddings
/// @details    Includes methods for computing membrane embeddings from search,
///    sequence, structure, user input, and smaller component embeddings
///    Last Modified: 7/24/14
///
/// @author  Julia Koehler (julia.koehler1982@gmail.com)
/// @author  Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_geometry_Embedding_fwd_hh
#define INCLUDED_protocols_membrane_geometry_Embedding_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {
namespace geometry {

class Embedding;
typedef utility::pointer::shared_ptr< Embedding > EmbeddingOP;
typedef utility::pointer::shared_ptr< Embedding const > EmbeddingCOP;

} // geometry
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_geometry_Embedding_fwd_hh
