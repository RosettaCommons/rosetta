// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphContainer.cc
/// @brief A container for the BuriedUnsatPenaltyGraph, to allow it to be cached in a pose while skirting multiple inheritance issues.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#include <core/pack/guidance_scoreterms/buried_unsat_penalty/graph/BuriedUnsatPenaltyGraphContainer.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "core.pack.guidance_scoreterms.buried_unsat_penalty.graph.BuriedUnsatPenaltyGraphContainer" );


namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace buried_unsat_penalty {
namespace graph {

/// @brief Constructor must be provided with an OP to a graph, which it stores.
BuriedUnsatPenaltyGraphContainer::BuriedUnsatPenaltyGraphContainer( BuriedUnsatPenaltyGraphOP graph ):
	basic::datacache::CacheableData(),
	graph_( graph )
{
	runtime_assert( graph != nullptr );
}

/// @brief Destructor
BuriedUnsatPenaltyGraphContainer::~BuriedUnsatPenaltyGraphContainer(){}

/// @brief Copy constructor shallow-copies the graphOP
BuriedUnsatPenaltyGraphContainer::BuriedUnsatPenaltyGraphContainer( BuriedUnsatPenaltyGraphContainer const &src ) :
	basic::datacache::CacheableData( src ),
	graph_( src.graph_ )
{}

/// @brief Copy this object and return an owning pointer to the copy.
basic::datacache::CacheableDataOP
BuriedUnsatPenaltyGraphContainer::clone() const {
	return basic::datacache::CacheableDataOP( BuriedUnsatPenaltyGraphContainerOP( new BuriedUnsatPenaltyGraphContainer( *this ) ) );
}

} //core
} //pack
} //guidance_scoreterms
} //buried_unsat_penalty
} //graph






