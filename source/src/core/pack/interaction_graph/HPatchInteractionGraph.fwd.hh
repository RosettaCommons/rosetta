// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/HPatchInteractionGraph.fwd.hh
/// @brief  Forward header for an IG that designs against hydrophobic area on the surface of proteins
/// @author Ron Jacak

#ifndef INCLUDED_core_pack_interaction_graph_HPatchInteractionGraph_fwd_hh
#define INCLUDED_core_pack_interaction_graph_HPatchInteractionGraph_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Package headers
#include <core/pack/interaction_graph/PDInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.fwd.hh>

namespace core {
namespace pack {
namespace interaction_graph {

template < typename V, typename E, typename G > class AdditionalBackgroundNodesInteractionGraph;
template < typename V, typename E, typename G > class HPatchInteractionGraph;

typedef HPatchInteractionGraph< PDNode, PDEdge, PDInteractionGraph > PDHPatchInteractionGraph;
typedef HPatchInteractionGraph< LinearMemNode, LinearMemEdge, LinearMemoryInteractionGraph > LinearMemoryHPatchInteractionGraph;

typedef utility::pointer::shared_ptr< PDHPatchInteractionGraph > PDHPatchInteractionGraphOP;
typedef utility::pointer::shared_ptr< LinearMemoryHPatchInteractionGraph > LinearMemoryHPatchInteractionGraphOP;


}
}
}

#endif
