// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/graph/Digraph.fwd.hh
/// @brief  directed graph base class forward declarations
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_graph_Digraph_fwd_hh
#define INCLUDED_utility_graph_Digraph_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace utility {
namespace graph {

class DirectedEdgeListElement;
class DirectedEdgeList;
class DirectedEdgeListIterator;
class DirectedEdgeListConstIterator;

class DirectedNode;
class DirectedEdge;
class Digraph;

typedef utility::pointer::shared_ptr< Digraph > DigraphOP;
typedef utility::pointer::shared_ptr< Digraph const > DigraphCOP;

} // namespace graph
} // namespace utility

#endif
