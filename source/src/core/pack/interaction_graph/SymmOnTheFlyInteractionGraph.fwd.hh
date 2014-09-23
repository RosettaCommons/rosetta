// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/SymmOnTheFlyInteractionGraph.fwd.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_interaction_graph_SymmOnTheFlyInteractionGraph_fwd_hh
#define INCLUDED_core_pack_interaction_graph_SymmOnTheFlyInteractionGraph_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace interaction_graph {

class SymmOnTheFlyNode;
class SymmOnTheFlyEdge;
class SymmOnTheFlyInteractionGraph;

typedef utility::pointer::shared_ptr< SymmOnTheFlyInteractionGraph > SymmOnTheFlyInteractionGraphOP;
typedef utility::pointer::shared_ptr< SymmOnTheFlyInteractionGraph const > SymmOnTheFlyInteractionGraphCOP;

}
}
}

#endif //INCLUDED_core_pack_interaction_graph_SymmOnTheFlyInteractionGraph_FWD_HH

