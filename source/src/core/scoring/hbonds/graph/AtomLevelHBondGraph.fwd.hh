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
/// @author Jack Maguire, jack@med.unc.edu


#ifndef INCLUDED_core_scoring_hbonds_graph_AtomLevelHBondGraph_FWD_HH
#define INCLUDED_core_scoring_hbonds_graph_AtomLevelHBondGraph_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <core/scoring/hbonds/graph/HBondGraph.fwd.hh>

namespace core {
namespace scoring {
namespace hbonds {
namespace graph {

class AtomLevelHBondNode;
typedef utility::pointer::shared_ptr< AtomLevelHBondNode > AtomLevelHBondNodeOP;
typedef utility::pointer::shared_ptr< AtomLevelHBondNode const > AtomLevelHBondNodeCOP;

class AtomLevelHBondEdge;
typedef utility::pointer::shared_ptr< AtomLevelHBondEdge > AtomLevelHBondEdgeOP;
typedef utility::pointer::shared_ptr< AtomLevelHBondEdge const > AtomLevelHBondEdgeCOP;

class AtomLevelHBondGraph;
typedef utility::pointer::shared_ptr< AtomLevelHBondGraph > AtomLevelHBondGraphOP;
typedef utility::pointer::shared_ptr< AtomLevelHBondGraph const > AtomLevelHBondGraphCOP;

} //graph
} //hbonds
} //scoring
} //core

#endif
