// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphContainer.hh
/// @brief A container for the BuriedUnsatPenaltyGraph, to allow it to be cached in a pose while skirting multiple inheritance issues.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


#ifndef INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_graph_BuriedUnsatPenaltyGraphContainer_hh
#define INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_graph_BuriedUnsatPenaltyGraphContainer_hh

#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphContainer.fwd.hh>
#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraph.fwd.hh>

#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace buried_unsat_penalty {
namespace graph {

///@brief A container for the BuriedUnsatPenaltyGraph, to allow it to be cached in a pose while skirting multiple inheritance issues.
class BuriedUnsatPenaltyGraphContainer : public basic::datacache::CacheableData {

public:

	/// @brief Default constructor is explicitly deleted.
	BuriedUnsatPenaltyGraphContainer() = delete;

	/// @brief Constructor must be provided with an OP to a graph, which it stores.
	BuriedUnsatPenaltyGraphContainer( BuriedUnsatPenaltyGraphOP graph );

	/// @brief Copy constructor shallow-copies the graphOP.
	BuriedUnsatPenaltyGraphContainer(BuriedUnsatPenaltyGraphContainer const & src);

	/// @brief Destructor
	virtual ~BuriedUnsatPenaltyGraphContainer();

	/// @brief Copy this object and return an owning pointer to the copy.
	basic::datacache::CacheableDataOP
	clone() const override;

	/// @brief Access the graph.
	inline BuriedUnsatPenaltyGraphOP graph() { return graph_; }

	/// @brief Access the graph (const-access).
	inline BuriedUnsatPenaltyGraphCOP graph() const { return graph_; }

private:

	/// @brief This container just contains a BuriedUnsatPenaltyGraphOP.
	BuriedUnsatPenaltyGraphOP graph_;

};


} //core
} //pack
} //guidance_scoreterms
} //buried_unsat_penalty
} //graph



#endif //INCLUDED_core_pack_guidance_scoreterms_buried_unsat_penalty_graph_BuriedUnsatPenaltyGraphContainer_hh
