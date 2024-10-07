// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file protocols/cyclic_peptide/crosslinker/Lanthionine_Helper.hh
/// @brief A crosslinker for lanthipeptides (lanthionine and methyllanthionine).
/// @author Clay Tydings (claiborne.w.tydings@vanderbilt.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_crosslinker_Lanthionine_Helper_HH
#define INCLUDED_protocols_cyclic_peptide_crosslinker_Lanthionine_Helper_HH

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/Lanthionine_Helper.fwd.hh>
#include <protocols/cyclic_peptide/crosslinker/CrosslinkerMoverHelper.hh>

// Protocol headers

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <utility/vector1.hh>

//--namespace--
namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

/// @brief A crosslinker class for lanthipeptides (lanthionine and methyllanthionine).
/// @author Clay Tydings (claiborne.w.tydings@vanderbilt.edu)
class Lanthionine_Helper : public protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelper {

public: //Constructors

	/// @brief Default constructor
	Lanthionine_Helper();

	/// @brief Copy constructor
	Lanthionine_Helper( Lanthionine_Helper const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~Lanthionine_Helper() override;

public: // public methods

	/// @brief Given a pose and a selection of residues, add the linker,
	/// align it crudely to the selected residues, and set up covalent bonds.
	/// @details Must be defined by derived classes.  Version for asymmetric poses.
	void
	add_linker_asymmetric(
		core::pose::Pose & pose,
		core::select::residue_selector::ResidueSubset const & selection
	) const override;

	/// @brief Given a pose and a linker, add bonds between the linker and the residues that
	/// coordinate the linker.
	/// @details Can be called by add_linker_asymmetric().  Must be defined by derived classes
	/// (pure virtual).  Version for asymmetric poses.
	void
	add_linker_bonds_asymmetric(
		core::pose::Pose & pose,
		utility::vector1< core::Size > const & res_indices,
		core::Size const linker_index
	) const override;

	/// @brief Given a pose and a selection of residues, add the linker,
	/// align it crudely to the selected residues, and set up covalent bonds.
	/// @details Must be defined by derived classes.  Version for symmetric poses.
	void
	add_linker_symmetric(
		core::pose::Pose & pose,
		core::select::residue_selector::ResidueSubset const & selection
	) const override;

	/// @brief Given a pose and a linker, add bonds between the linker and the residues
	/// that coordinate the linker.
	/// @details Can be called by add_linker_symmetric().  Must be defined by derived
	/// classes (pure virtual).  Version for symmetric poses.
	void
	add_linker_bonds_symmetric(
		core::pose::Pose & pose,
		core::Size const res1,
		core::Size const linker_index1,
		core::Size const linker_index2
	) const override;

	/// @brief Given a selection of residues that have already been connected to a crosslinker,
	/// add constraints for the crosslinker.
	/// @details Must be defined by derived classes.  Version for asymmetric poses.
	void
	add_linker_constraints_asymmetric(
		core::pose::Pose & pose,
		core::select::residue_selector::ResidueSubset const & selection
	) const override;

	/// @brief Given a selection of residues that have already been connected to a crosslinker,
	/// add constraints for the crosslinker.
	/// @details Must be defined by derived classes.  Version for symmetric poses.
	void
	add_linker_constraints_symmetric(
		core::pose::Pose & pose,
		core::select::residue_selector::ResidueSubset const & selection,
		bool const linker_was_added
	) const override;

	/// @brief Given indices of residues that are already linked to a linker, get the index
	/// of the linker.
	/// @details Throws an error if the residues are not all linked to the same linker
	/// residue.  Must be defined by derived classes.
	core::Size
	get_linker_index_asymmetric(
		core::pose::Pose const & pose,
		utility::vector1< core::Size > const & res_indices
	) const override;

	/// @brief Given indices of residues that are already linked to pieces of a linker, get
	/// of the indices of the symmetric pieces of the linker.
	/// @details Throws an error if a residue is not linked to something.  Must be defined by
	/// derived classes.
	void
	get_linker_indices_symmetric(
		core::pose::Pose const & pose,
		utility::vector1< core::Size > const & res_indices,
		utility::vector1< core::Size > & linker_indices
	) const override;

	/// @brief Given a pose with residues selected to be linked by a linker, determine whether
	/// the residues are too far apart.
	/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	bool
	filter_by_sidechain_distance_asymmetric(
		core::pose::Pose const & pose,
		core::select::residue_selector::ResidueSubset const & selection,
		core::Real const & filter_multiplier
	) const override;

	/// @brief Given a pose with residues selected to be linked by a linker, determine whether
	/// the residues are too far apart.
	/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
	/// This version is for symmetric poses.
	/// @note Higher values of the filter multiplier make it more permissive.
	bool
	filter_by_sidechain_distance_symmetric(
		core::pose::Pose const & pose,
		core::select::residue_selector::ResidueSubset const & selection,
		core::Real const &filter_multiplier
	) const override;

	/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
	/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	bool
	filter_by_constraints_energy_asymmetric(
		core::pose::Pose const & pose,
		core::select::residue_selector::ResidueSubset const & selection,
		core::Real const & filter_multiplier
	) const override;

	/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
	/// This version is for symmetric poses.
	/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	bool
	filter_by_constraints_energy_symmetric(
		core::pose::Pose const & pose,
		core::select::residue_selector::ResidueSubset const & selection,
		bool const linker_was_added, core::Real const & filter_multiplier
	) const override;

	/// @brief Does this CrosslinkerMoverHelper add a residue for the linker?
	bool
	helper_adds_linker_residue() const override {
		return false;
	}

private: // private methods

	///@brief Given a selection of exactly two residues and a pose, return the two lanthionine
	///residues.

	void
	get_lanthionine_residues_from_selection(
		core::select::residue_selector::ResidueSubset const & selection,
		core::pose::Pose const & pose,
		core::Size & dalares,
		core::Size & cysres
	) const;


private: // data

};

//--end_namespace--
} //crosslinker
} //cyclic_peptide
} //protocols

#endif //protocols_cyclic_peptide_crosslinker_Lanthionine_Helper_HH
