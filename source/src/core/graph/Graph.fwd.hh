// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/graph/Graph.fwd.hh
/// @brief  graph base class forward declarations
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_graph_Graph_fwd_hh
#define INCLUDED_core_graph_Graph_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace graph {

class EdgeListElement;
class EdgeList;
class EdgeListIterator;
class EdgeListConstIterator;

class Node;
class Edge;
class Graph;

typedef utility::pointer::shared_ptr< Graph > GraphOP;
typedef utility::pointer::shared_ptr< Graph const > GraphCOP;

} // namespace graph
} // namespace core

#endif
