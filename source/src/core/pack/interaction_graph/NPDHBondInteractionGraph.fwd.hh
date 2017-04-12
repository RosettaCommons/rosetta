// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/NPDHBondInteractionGraph.fwd.hh
/// @brief  Forward header for an IG that assigns hydrogen bonds a (non-pairwise-decomposable) weight
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_interaction_graph_NPDHBondInteractionGraph_fwd_hh
#define INCLUDED_core_pack_interaction_graph_NPDHBondInteractionGraph_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Package headers
#include <core/pack/interaction_graph/PDInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/DensePDInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.fwd.hh>

namespace core {
namespace pack {
namespace interaction_graph {

struct NPDHBond;
typedef utility::pointer::shared_ptr< NPDHBond > NPDHBondOP;

template < typename V, typename E, typename G > class AdditionalBackgroundNodesInteractionGraph;
template < typename V, typename E, typename G > class NPDHBondInteractionGraph;

// For use in design, but perhaps, slower generally than the Linmem version
typedef NPDHBondInteractionGraph< PDNode, PDEdge, PDInteractionGraph > StandardNPDHBondInteractionGraph;
// For use in fixed-sequence repacking
typedef NPDHBondInteractionGraph< DensePDNode, DensePDEdge, DensePDInteractionGraph > DenseNPDHBondInteractionGraph;
// For use in design
typedef NPDHBondInteractionGraph< LinearMemNode, LinearMemEdge, LinearMemoryInteractionGraph > LinearMemoryNPDHBondInteractionGraph;

typedef utility::pointer::shared_ptr< StandardNPDHBondInteractionGraph > StandardNPDHBondInteractionGraphOP;
typedef utility::pointer::shared_ptr< DenseNPDHBondInteractionGraph > DenseNPDHBondInteractionGraphOP;
typedef utility::pointer::shared_ptr< LinearMemoryNPDHBondInteractionGraph > LinearMemoryNPDHBondInteractionGraphOP;


}
}
}

#endif
