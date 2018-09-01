// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file utility/graph/LowMemGraph.fwd.hh
/// @brief A lower memory version of utility::graph::Graph with three key limitations
///        1. Due to std::vector::resize(), all of the Edge* can go invalid any time you call add_edge().
///        2. This doesn't have constant time deletion (a definite design goal for utility::graph::Graph)
///        3. Deleting an element doesn't actually delete it from memory
/// @author Brian Coventry (bcov@uw.edu)

#ifndef INCLUDED_utility_graph_LowMemGraph_FWD_HH
#define INCLUDED_utility_graph_LowMemGraph_FWD_HH

#include <utility/pointer/owning_ptr.hh>


namespace utility {
namespace graph {



class LowMemNode;
class LowMemEdge;
class LowMemEdgeListIter;
class LowMemEdgeListConstIter;

// Please don't define LowMemNodeOP and LowMemEdgeOP
// These classes should have virtual destructors but don't.
// In the normal use of this class, one should never have anything
// besides LowMemNode* and LowMemEdge*, so this should
// never matter. Defining the OPs will make it easier for people
// to make a mistake.

template< class _Node, class _Edge >
class LowMemGraph;

class LowMemGraphBase;

typedef LowMemGraph< LowMemNode, LowMemEdge > DefaultLowMemGraph;

typedef utility::pointer::shared_ptr< DefaultLowMemGraph > DefaultLowMemGraphOP;
typedef utility::pointer::shared_ptr< DefaultLowMemGraph const > DefaultLowMemGraphCOP;



}
}

#endif
