// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/HierarchicalHybridJD_ResultsSummaryBase.hh
/// @brief A pure virtual base class for the helper classes used by the HierarchicalHybridJDApplication
/// class for transmitting information up the MPI hierarchy.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_ResultsSummaryBase_hh
#define INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_ResultsSummaryBase_hh

#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_ResultsSummaryBase.fwd.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace cyclic_peptide_predict {

/// @brief A pure virtual base class for the helper classes used by the HierarchicalHybridJDApplication
/// class for transmitting information up the MPI hierarchy.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class HierarchicalHybridJD_ResultsSummaryBase : public utility::VirtualBase {

public:

	/// @brief Default constructor.
	HierarchicalHybridJD_ResultsSummaryBase() = default;

	/// @brief Copy constructor.
	HierarchicalHybridJD_ResultsSummaryBase(HierarchicalHybridJD_ResultsSummaryBase const & ) = default;

	/// @brief Constructor with options
	///
	HierarchicalHybridJD_ResultsSummaryBase(
		int const originating_node_MPI_rank,
		core::Size const jobindex_on_originating_node
	);

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	virtual HierarchicalHybridJD_ResultsSummaryBaseOP clone() const = 0;

public:
	/// ------------- Getters -------------------------------

	/// @brief Get the MPI rank of the process that initially carried out this job.
	///
	inline int originating_node_MPI_rank() const { return originating_node_MPI_rank_; }

	/// @brief Get the local index of the job on the node on which it was carried out.
	///
	inline core::Size jobindex_on_originating_node() const { return jobindex_on_originating_node_; }

	/// @brief Get the list of procs that transmitted this JobResultsSummary.
	/// @details Whenever a proc sends it to someone else, the receiving proc appends the sender's rank here.  This allows the director
	/// node to address an output request to the correct proc.
	inline utility::vector1 <int> const & MPI_ranks_handling_message() const { return mpi_ranks_handling_message_; }

public:
	/// ------------- Setters -------------------------------

	/// @brief Set the MPI rank of the process that initially carried out this job.
	///
	void set_originating_node_MPI_rank( int const rank);

	/// @brief Set the local index of the job on the node on which it was carried out.
	///
	void set_jobindex_on_originating_node( core::Size const index);

	/// @brief Add a proc that handled this message to the list of procs that transmitted this JobResultsSummary.
	/// @details Whenever a proc sends it to someone else, the receiving proc appends the sender's rank here.  This allows the director
	/// node to address an output request to the correct proc.
	void add_MPI_rank_handling_message( int const new_proc );

private:

	/// @brief The rank of the process that actually ran the job.
	///
	int originating_node_MPI_rank_ = 0;

	/// @brief The local index of the job on the process that ran the job.
	///
	core::Size jobindex_on_originating_node_ = 0;

	/// @brief The list of procs that transmitted this JobResultsSummary.
	/// @details Whenever a proc sends it to someone else, the receiving proc appends the sender's rank here.  This allows the director
	/// node to address an output request to the correct proc.
	utility::vector1 < int > mpi_ranks_handling_message_;

};

} //cyclic_peptide_predict
} //protocols

#endif //INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_ResultsSummaryBase_hh
