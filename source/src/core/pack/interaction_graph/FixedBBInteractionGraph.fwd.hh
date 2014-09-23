// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/FixedBBInteractionGraph.fwd.hh
/// @brief  Interaction graph base class for fixed-backbone packing
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_interaction_graph_FixedBBInteractionGraph_fwd_hh
#define INCLUDED_core_pack_interaction_graph_FixedBBInteractionGraph_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace interaction_graph {

class FixedBBNode;
class FixedBBEdge;
class FixedBBInteractionGraph;

typedef utility::pointer::shared_ptr< FixedBBInteractionGraph > FixedBBInteractionGraphOP;
typedef utility::pointer::shared_ptr< FixedBBInteractionGraph const > FixedBBInteractionGraphCOP;

} //end namespace interaction_graph
} //end namespace pack
} //end namespace core

#endif
