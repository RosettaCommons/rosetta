// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/HierarchicalHybridJD_PNearToArbitraryStateSummary.hh
/// @brief A class for storing the PNear, Keq, and DeltaG_folding values to an arbitrary state that has been sampled.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_PNearToArbitraryStateSummary_hh
#define INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_PNearToArbitraryStateSummary_hh

#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_PNearToArbitraryStateSummary.fwd.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace protocols {
namespace cyclic_peptide_predict {

/// @brief A class for storing the PNear, Keq, and DeltaG_folding values to an arbitrary state that has been sampled.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class HierarchicalHybridJD_PNearToArbitraryStateSummary : public utility::VirtualBase {

public:

	/// @brief Default constructor (explicitly deleted).
	HierarchicalHybridJD_PNearToArbitraryStateSummary() = delete;

	/// @brief Options constructor.
	HierarchicalHybridJD_PNearToArbitraryStateSummary(
		core::Real const pnear,
		core::Real const Keq,
		core::Real const deltaG_folding,
		core::Size const originating_mpi_node,
		core::Size const jobindex_on_originating_node
	);

	/// @brief Destructor.
	~HierarchicalHybridJD_PNearToArbitraryStateSummary() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	HierarchicalHybridJD_PNearToArbitraryStateSummaryOP clone() const;

public: //Accessors

	/// @brief Access the PNear value.
	inline core::Real pnear() const { return pnear_; }

	/// @brief Access the PNfolding equilibtrium constant value.
	inline core::Real Keq() const { return Keq_; }

	/// @brief Access the Delta-G of folding value.
	inline core::Real deltaG_folding() const { return deltaG_folding_; }

	/// @brief Access the index of the originating MPI node.
	inline core::Size originating_mpi_node() const { return originating_mpi_node_; }

	/// @brief Access the job index on the originating node.
	inline core::Size jobindex_on_originating_node() const { return jobindex_on_originating_node_; }

private:

	/// @brief The PNear value.
	core::Real pnear_ = 0.0;

	/// @brief The folding equilbrium constant.
	core::Real Keq_ = 0.0;

	/// @brief The Delta-G of folding.
	core::Real deltaG_folding_ = 0.0;

	/// @brief The MPI node that originated this structure.
	core::Size originating_mpi_node_ = 0;

	/// @brief The job index on the originating node.
	core::Size jobindex_on_originating_node_ = 0;

};

} //cyclic_peptide_predict
} //protocols

#endif //INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_PNearToArbitraryStateSummary_hh
