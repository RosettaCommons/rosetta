// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/membrane/geometry/EmbeddingDef.fwd.hh
///
/// @brief      Basic Embedding Definitions for Membrane Embedding
/// @details    Class contains a normal and center
///    Last Modified: 6/17/14
///
/// @author  JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_geometry_EmbeddingDef_fwd_hh
#define INCLUDED_protocols_membrane_geometry_EmbeddingDef_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {
namespace geometry {

class EmbeddingDef;
typedef utility::pointer::shared_ptr< EmbeddingDef > EmbeddingDefOP;
typedef utility::pointer::shared_ptr< EmbeddingDef const > EmbeddingDefCOP;

} // geometry
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_geometry_EmbeddingDef_fwd_hh


