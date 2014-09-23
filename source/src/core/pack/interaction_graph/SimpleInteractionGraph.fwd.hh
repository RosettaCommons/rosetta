// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/SimpleInteractionGraph.fwd.hh
/// @brief
/// @author Liz Kellogg (ekellogg@u.washington.edu)

#ifndef INCLUDED_core_pack_interaction_graph_SimpleInteractionGraph_FWD_HH
#define INCLUDED_core_pack_interaction_graph_SimpleInteractionGraph_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace interaction_graph {

class SimpleNode;
class SimpleEdge;
class SimpleInteractionGraph;

typedef utility::pointer::shared_ptr< SimpleInteractionGraph > SimpleInteractionGraphOP;
typedef utility::pointer::shared_ptr< SimpleInteractionGraph const > SimpleInteractionGraphCOP;

}
}
}

#endif


