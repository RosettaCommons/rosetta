// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/hbonds/graph/HBondGraph.fwd.hh
/// @brief AtomLevelHBondGraph forward declarations
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_core_scoring_hbonds_graph_AtomLevelHBondGraph_FWD_HH
#define INCLUDED_core_scoring_hbonds_graph_AtomLevelHBondGraph_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace hbonds {
namespace graph {

class AtomLevelHBondNode;
using AtomLevelHBondNodeOP = utility::pointer::shared_ptr< AtomLevelHBondNode >;
using AtomLevelHBondNodeCOP = utility::pointer::shared_ptr< AtomLevelHBondNode const >;

class AtomLevelHBondEdge;
using AtomLevelHBondEdgeOP = utility::pointer::shared_ptr< AtomLevelHBondEdge >;
using AtomLevelHBondEdgeCOP = utility::pointer::shared_ptr< AtomLevelHBondEdge const >;

class AtomLevelHBondGraph;
using AtomLevelHBondGraphOP = utility::pointer::shared_ptr< AtomLevelHBondGraph >;
using AtomLevelHBondGraphCOP = utility::pointer::shared_ptr< AtomLevelHBondGraph const >;

} //graph
} //hbonds
} //scoring
} //core

#endif
