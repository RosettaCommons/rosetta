// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/HierarchicalHybridJD_RMSDToBestSummary.hh
/// @brief A small helper class used by the HierarchicalHybridJDApplication class for transmitting RMSD information up the MPI hierarchy.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_RMSDToBestSummary_hh
#define INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_RMSDToBestSummary_hh

#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_RMSDToBestSummary.fwd.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_ResultsSummaryBase.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace protocols {
namespace cyclic_peptide_predict {

/// @brief A small helper class used by the HierarchicalHybridJDApplication class for transmitting RMSD information up the MPI hierarchy.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class HierarchicalHybridJD_RMSDToBestSummary : public HierarchicalHybridJD_ResultsSummaryBase {

public:

	/// @brief Default constructor.
	HierarchicalHybridJD_RMSDToBestSummary() = default;

	/// @brief Options constructor/
	HierarchicalHybridJD_RMSDToBestSummary(
		int const originating_node_MPI_rank,
		core::Size const jobindex_on_originating_node,
		core::Real const & rmsd_to_best
	);

	/// @brief Copy constructor.
	HierarchicalHybridJD_RMSDToBestSummary(HierarchicalHybridJD_RMSDToBestSummary const & ) = default;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	HierarchicalHybridJD_ResultsSummaryBaseOP clone() const override;

public:
	//Getters

	/// @brief Get the RMSD to the best structure.
	inline core::Real const & rmsd_to_best() const { return rmsd_to_best_; }

public:
	//Setters

	inline void set_rmsd_to_best( core::Real const & setting ) { rmsd_to_best_ = setting; }

private:

	/// @brief The RMSD to the best structure.
	core::Real rmsd_to_best_ = 0.0;

};

} //cyclic_peptide_predict
} //protocols

#endif //INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_RMSDToBestSummary_hh
