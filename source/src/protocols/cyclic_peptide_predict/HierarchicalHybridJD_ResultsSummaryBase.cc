// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/HierarchicalHybridJD_ResultsSummaryBase.cc
/// @brief A pure virtual base class for the helper classes used by the HierarchicalHybridJDApplication
/// class for transmitting information up the MPI hierarchy.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_ResultsSummaryBase.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.cyclic_peptide_predict.HierarchicalHybridJD_ResultsSummaryBase" );


namespace protocols {
namespace cyclic_peptide_predict {


/// @brief Constructor with options
///
HierarchicalHybridJD_ResultsSummaryBase::HierarchicalHybridJD_ResultsSummaryBase(
	int const originating_node_MPI_rank,
	core::Size const jobindex_on_originating_node
) :
	originating_node_MPI_rank_( originating_node_MPI_rank ),
	jobindex_on_originating_node_( jobindex_on_originating_node )
{
	runtime_assert( originating_node_MPI_rank_ > 0 );
	runtime_assert( jobindex_on_originating_node_ > 0);
}

/// @brief Set the MPI rank of the process that initially carried out this job.
///
void
HierarchicalHybridJD_ResultsSummaryBase::set_originating_node_MPI_rank(
	int const rank
) {
	originating_node_MPI_rank_ = rank;
}

/// @brief Set the local index of the job on the node on which it was carried out.
///
void
HierarchicalHybridJD_ResultsSummaryBase::set_jobindex_on_originating_node(
	core::Size const index
) {
	jobindex_on_originating_node_ = index;
}

/// @brief Add a proc that handled this message to the list of procs that transmitted this JobResultsSummary.
/// @details Whenever a proc sends it to someone else, the receiving proc appends the sender's rank here.  This allows the director
/// node to address an output request to the correct proc.
void
HierarchicalHybridJD_ResultsSummaryBase::add_MPI_rank_handling_message(
	int const new_proc
) {
	mpi_ranks_handling_message_.push_back( new_proc );
}

} //cyclic_peptide_predict
} //protocols
