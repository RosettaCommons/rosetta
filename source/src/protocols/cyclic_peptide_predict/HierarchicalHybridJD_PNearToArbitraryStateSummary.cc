// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/HierarchicalHybridJD_PNearToArbitraryStateSummary.cc
/// @brief A class for storing the PNear, Keq, and DeltaG_folding values to an arbitrary state that has been sampled.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_PNearToArbitraryStateSummary.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.cyclic_peptide_predict.HierarchicalHybridJD_PNearToArbitraryStateSummary" );


namespace protocols {
namespace cyclic_peptide_predict {

/// @brief Options constructor.
HierarchicalHybridJD_PNearToArbitraryStateSummary::HierarchicalHybridJD_PNearToArbitraryStateSummary(
	core::Real const pnear,
	core::Real const Keq,
	core::Real const deltaG_folding,
	core::Size const originating_mpi_node,
	core::Size const jobindex_on_originating_node
) {
	pnear_ = pnear;
	Keq_ = Keq;
	deltaG_folding_ = deltaG_folding;
	originating_mpi_node_ = originating_mpi_node;
	jobindex_on_originating_node_ = jobindex_on_originating_node;
}

/// @brief Destructor.
HierarchicalHybridJD_PNearToArbitraryStateSummary::~HierarchicalHybridJD_PNearToArbitraryStateSummary() = default;

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
HierarchicalHybridJD_PNearToArbitraryStateSummaryOP
HierarchicalHybridJD_PNearToArbitraryStateSummary::clone() const {
	return utility::pointer::make_shared< HierarchicalHybridJD_PNearToArbitraryStateSummary >( *this );
}

} //cyclic_peptide_predict
} //protocols
