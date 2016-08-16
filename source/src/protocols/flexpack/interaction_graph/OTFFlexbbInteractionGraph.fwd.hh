// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/flexpack/
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_flexpack_interaction_graph_OTFFlexbbInteractionGraph_fwd_hh
#define INCLUDED_protocols_flexpack_interaction_graph_OTFFlexbbInteractionGraph_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace flexpack {
namespace interaction_graph {

class OTFFlexbbNode;
class OTFFlexbbEdge;
class OTFFlexbbInteractionGraph;

typedef utility::pointer::shared_ptr< OTFFlexbbInteractionGraph > OTFFlexbbInteractionGraphOP;
typedef utility::pointer::shared_ptr< OTFFlexbbInteractionGraph const > OTFFlexbbInteractionGraphCOP;

}
}
}

#endif


