// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/conformation/AtomGraph.fwd.hh
/// @author Sam DeLuca

#ifndef INCLUDED_core_conformation_AtomGraph_fwd_hh
#define INCLUDED_core_conformation_AtomGraph_fwd_hh

#include <core/conformation/AtomGraphData.fwd.hh>
#include <utility/graph/UpperEdgeGraph.fwd.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace conformation {

typedef utility::graph::UpperEdgeGraph<AtomGraphVertexData, AtomGraphEdgeData> AtomGraph;
typedef utility::pointer::shared_ptr<AtomGraph > AtomGraphOP;

}
}

#endif /* ATOMGRAPH_FWD_HH_ */
