// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/hbonds/graph/HBondGraph.fwd.hh
/// @brief HBondGraph (derived from core::graph::Graph) forward declarations
/// @author Jack Maguire, jackmaguire1444@gmail.com


#ifndef INCLUDED_core_scoring_hbonds_graph_HBondGraph_FWD_HH
#define INCLUDED_core_scoring_hbonds_graph_HBondGraph_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace hbonds {
namespace graph {

class HBondNode;
typedef utility::pointer::shared_ptr< HBondNode > HBondNodeOP;
typedef utility::pointer::shared_ptr< HBondNode const > HBondNodeCOP;

class HBondEdge;
typedef utility::pointer::shared_ptr< HBondEdge > HBondEdgeOP;
typedef utility::pointer::shared_ptr< HBondEdge const > HBondEdgeCOP;

class AbstractHBondGraph;
typedef utility::pointer::shared_ptr< AbstractHBondGraph > AbstractHBondGraphOP;
typedef utility::pointer::shared_ptr< AbstractHBondGraph const > AbstractHBondGraphCOP;

class HBondGraph;
typedef utility::pointer::shared_ptr< HBondGraph > HBondGraphOP;
typedef utility::pointer::shared_ptr< HBondGraph const > HBondGraphCOP;

} //graph
} //hbonds
} //scoring
} //core

#endif
