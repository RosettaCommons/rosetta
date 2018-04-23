// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphContainer.fwd.hh
/// @brief A container for the BuriedUnsatPenaltyGraph, to allow it to be cached in a pose while skirting multiple inheritance issues.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_graph_BuriedUnsatPenaltyGraphContainer_fwd_hh
#define INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_graph_BuriedUnsatPenaltyGraphContainer_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace buried_unsat_penalty {
namespace graph {

class BuriedUnsatPenaltyGraphContainer;

typedef utility::pointer::shared_ptr< BuriedUnsatPenaltyGraphContainer > BuriedUnsatPenaltyGraphContainerOP;
typedef utility::pointer::shared_ptr< BuriedUnsatPenaltyGraphContainer const > BuriedUnsatPenaltyGraphContainerCOP;

} //core
} //pack
} //guidance_scoreterms
} //buried_unsat_penalty
} //graph

#endif //INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_graph_BuriedUnsatPenaltyGraphContainer_fwd_hh
