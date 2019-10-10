// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/HierarchicalHybridJD_RMSDToBestSummary.cc
/// @brief A small helper class used by the HierarchicalHybridJDApplication class for transmitting RMSD information up the MPI hierarchy.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_RMSDToBestSummary.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.cyclic_peptide_predict.HierarchicalHybridJD_RMSDToBestSummary" );


namespace protocols {
namespace cyclic_peptide_predict {

/// @brief Options constructor/
HierarchicalHybridJD_RMSDToBestSummary::HierarchicalHybridJD_RMSDToBestSummary(
	int const originating_node_MPI_rank,
	core::Size const jobindex_on_originating_node,
	core::Real const & rmsd_to_best
) :
	HierarchicalHybridJD_ResultsSummaryBase( originating_node_MPI_rank, jobindex_on_originating_node ),
	rmsd_to_best_( rmsd_to_best )
{}

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
HierarchicalHybridJD_ResultsSummaryBaseOP
HierarchicalHybridJD_RMSDToBestSummary::clone() const {
	return utility::pointer::make_shared< HierarchicalHybridJD_RMSDToBestSummary >( *this );
}

} //cyclic_peptide_predict
} //protocols
