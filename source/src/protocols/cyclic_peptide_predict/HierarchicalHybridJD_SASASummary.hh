// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/HierarchicalHybridJD_SASASummary.hh
/// @brief A class for storing information about solvent-exposed surface area and for transmitting it up the MPI hierarchy.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_SASASummary_hh
#define INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_SASASummary_hh

#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_SASASummary.fwd.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_ResultsSummaryBase.hh>

// Corew headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace protocols {
namespace cyclic_peptide_predict {

/// @brief A class for storing information about solvent-exposed surface area and for transmitting it up the MPI hierarchy.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class HierarchicalHybridJD_SASASummary : public HierarchicalHybridJD_ResultsSummaryBase {

public:

	/// @brief Default constructor.
	HierarchicalHybridJD_SASASummary() = default;

	/// @brief Options constructor.
	HierarchicalHybridJD_SASASummary(
		int const originating_node_MPI_rank,
		core::Size const jobindex_on_originating_node,
		core::Real const &sasa,
		core::Real const &polar_sasa,
		core::Real const &hydrophobic_sasa
	);

	/// @brief Copy constructor.
	HierarchicalHybridJD_SASASummary(HierarchicalHybridJD_SASASummary const & ) = default;

	/// @brief Destructor.
	~HierarchicalHybridJD_SASASummary() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	HierarchicalHybridJD_ResultsSummaryBaseOP clone() const override;

public: //Setters

	/// @brief Set total sasa, polar sasa, and hydrophobic sasa.
	/// @details Computes fraction_polar_sasa_ and fraction_hydrophobic_sasa_.
	void set(
		core::Real const &sasa,
		core::Real const &polar_sasa,
		core::Real const &hydrophobic_sasa
	);

public: //Getters

	/// @brief Get the total solvent-accessible surface area, in square Angstroms.
	inline core::Real const & sasa() const { return sasa_; }

	/// @brief Get the polar solvent-accessible surface area, in square Angstroms.
	inline core::Real const & polar_sasa() const { return polar_sasa_; }

	/// @brief Get the hydrophobic solvent-accessible surface area, in square Angstroms.
	inline core::Real const & hydrophobic_sasa() const { return hydrophobic_sasa_; }

	/// @brief Get the solvent-accessible polar surface area, as a fraction of the total surface area.
	inline core::Real const & fraction_polar_sasa() const { return fraction_polar_sasa_; }

	/// @brief Get the solvent-accessible hydrophobic surface area, as a fraction of the total surface area.
	inline core::Real const & fraction_hydrophobic_sasa() const { return fraction_hydrophobic_sasa_; }

private:

	/// @brief The total solvent-accessible surface area, in square Angstroms.
	core::Real sasa_ = 0.0;

	/// @brief The total polar solvent-accessible surface area, in square Angstroms.
	core::Real polar_sasa_ = 0.0;

	/// @brief The total hydrophobic solvent-accessible surface area, in square Angstroms.
	core::Real hydrophobic_sasa_ = 0.0;

	/// @brief The solvent-accessible polar surface area, as a fraction of the total surface area.
	core::Real fraction_polar_sasa_ = 0.0;

	/// @brief The solvent-accessible hydrophobic surface area, as a fraction of the total surface area.
	core::Real fraction_hydrophobic_sasa_ = 0.0;

};

} //cyclic_peptide_predict
} //protocols

#endif //INCLUDED_protocols_cyclic_peptide_predict_HierarchicalHybridJD_SASASummary_hh
