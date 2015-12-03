// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief      Optimizes the protein embedding in the membrane
/// @details Optimizes the protein embedding in the membrane given the smooth
///   high-res score function; transforms the protein into the membrane,
///   optimizes the membrane position (flexible), and uses the optimized
///   embedding to reposition the protein in the membrane
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_OptimizeProteinEmbeddingMover_fwd_hh
#define INCLUDED_protocols_membrane_OptimizeProteinEmbeddingMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace membrane {

class OptimizeProteinEmbeddingMover;
typedef utility::pointer::shared_ptr< OptimizeProteinEmbeddingMover > OptimizeProteinEmbeddingMoverOP;
typedef utility::pointer::shared_ptr< OptimizeProteinEmbeddingMover const > OptimizeProteinEmbeddingMoverCOP;

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_OptimizeProteinEmbeddingMover_fwd_hh
