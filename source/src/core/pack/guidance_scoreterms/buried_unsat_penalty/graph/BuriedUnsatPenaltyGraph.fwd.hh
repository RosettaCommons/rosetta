// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraph.fwd.hh
/// @brief Forward declarations for the BuriedUnsatPenaltyGraph class and its related Node and Edge classes.
/// @author Vikram K. Mulligan (vmullig@uw.edu).

#ifndef INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_graph_BuriedUnsatPenaltyGraph_FWD_HH
#define INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_graph_BuriedUnsatPenaltyGraph_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace buried_unsat_penalty {
namespace graph {

class BuriedUnsatPenaltyGraphHbond;
typedef utility::pointer::shared_ptr< BuriedUnsatPenaltyGraphHbond > BuriedUnsatPenaltyGraphHbondOP;
typedef utility::pointer::shared_ptr< BuriedUnsatPenaltyGraphHbond const > BuriedUnsatPenaltyGraphHbondCOP;

class BuriedUnsatPenaltyNodeData;
typedef utility::pointer::shared_ptr< BuriedUnsatPenaltyNodeData > BuriedUnsatPenaltyNodeDataOP;
typedef utility::pointer::shared_ptr< BuriedUnsatPenaltyNodeData const > BuriedUnsatPenaltyNodeDataCOP;

class BuriedUnsatPenaltyNode;
typedef utility::pointer::shared_ptr< BuriedUnsatPenaltyNode > BuriedUnsatPenaltyNodeOP;
typedef utility::pointer::shared_ptr< BuriedUnsatPenaltyNode const > BuriedUnsatPenaltyNodeCOP;

class BuriedUnsatPenaltyEdge;
typedef utility::pointer::shared_ptr< BuriedUnsatPenaltyEdge > BuriedUnsatPenaltyEdgeOP;
typedef utility::pointer::shared_ptr< BuriedUnsatPenaltyEdge const > BuriedUnsatPenaltyEdgeCOP;

class BuriedUnsatPenaltyGraph;
typedef utility::pointer::shared_ptr< BuriedUnsatPenaltyGraph > BuriedUnsatPenaltyGraphOP;
typedef utility::pointer::shared_ptr< BuriedUnsatPenaltyGraph const > BuriedUnsatPenaltyGraphCOP;

} //graph
} //buried_unsat_penalty
} //guidance terms
} //pack
} //core

#endif
