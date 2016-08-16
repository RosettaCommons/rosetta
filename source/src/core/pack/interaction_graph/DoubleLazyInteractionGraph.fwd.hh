// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/DoubleLazyInteractionGraph.fwd.hh
/// @brief  Interaction graph that computes each rotamer pair energy at most once
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_interaction_graph_DoubleLazyInteractionGraph_fwd_hh
#define INCLUDED_core_pack_interaction_graph_DoubleLazyInteractionGraph_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace interaction_graph {

class DoubleLazyNode;
class DoubleLazyEdge;
class DoubleLazyInteractionGraph;

typedef utility::pointer::shared_ptr< DoubleLazyInteractionGraph > DoubleLazyInteractionGraphOP;
typedef utility::pointer::shared_ptr< DoubleLazyInteractionGraph const > DoubleLazyInteractionGraphCOP;

}
}
}

#endif


