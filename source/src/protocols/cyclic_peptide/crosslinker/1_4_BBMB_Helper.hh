// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/1_4_BBMB_Helper.hh
/// @brief A derived class of the CrosslinkerMoverHelper base class, used to set up
/// the 1,4-bis(bromomethyl)benzene (BBMB) cross-linker.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_crosslinker_1_4_BBMB_Helper_hh
#define INCLUDED_protocols_cyclic_peptide_crosslinker_1_4_BBMB_Helper_hh

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/1_4_BBMB_Helper.fwd.hh>
#include <protocols/cyclic_peptide/crosslinker/CrosslinkerMoverHelper.hh>

// Protocol headers

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <core/types.hh>

namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

/// @brief A derived class of the CrosslinkerMoverHelper base class, used to set up
/// the 1,4-bis(bromomethyl)benzene (BBMB) cross-linker.
class One_Four_BBMB_Helper : public CrosslinkerMoverHelper {

public: //Constructors

	/// @brief Default constructor
	One_Four_BBMB_Helper();

	/// @brief Copy constructor
	One_Four_BBMB_Helper( One_Four_BBMB_Helper const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~One_Four_BBMB_Helper() override;


public: // public methods

	/// @brief Given a pose and a selection of exactly two residues, add the BBMB linker,
	/// align it crudely to the selected residues, and set up covalent bonds.
	void add_linker_asymmetric(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const override;

	/// @brief Given a pose and a linker, add bonds between the linker and the residues that coordinate the linker.
	/// @details Called by add_linker_asymmetric().  Version for asymmetric poses.
	void add_linker_bonds_asymmetric(core::pose::Pose &pose, utility::vector1< core::Size > const &res_indices, core::Size const linker_index ) const override;

	/// @brief Given a pose and a selection of exactly two residues, add the BBMB linker,
	/// align it crudely to the selected residues, and set up covalent bonds.
	/// @details Version for symmetric poses.
	void add_linker_symmetric(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const override;

	/// @brief Given a pose and a linker, add bonds between the BBMB linker and the residues that coordinate the linker.
	/// @details Called by add_linker_symmetric().  Version for symmetric poses.
	void add_linker_bonds_symmetric(core::pose::Pose &pose, core::Size const res1, core::Size const linker_index1, core::Size const linker_index2 ) const override;

	/// @brief Given a selection of exactly two residues that have already been connected to a 1,4-bis(bromomethyl)benzene crosslinker,
	/// add constraints for the crosslinker.
	void add_linker_constraints_asymmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const override;

	/// @brief Given a selection of exactly two residues that have already been connected to a 1,4-bis(bromomethyl)benzene crosslinker,
	/// add constraints for the crosslinker.  This version is for symmetric poses.
	void add_linker_constraints_symmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection, bool const linker_was_added) const override;

	/// @brief Given indices of two cysteine residues that are already linked to a BBMB, get the index
	/// of the BBMB residue.
	/// @details Throws an error if the two cysteines are not all linked to the same BBMB residue.
	core::Size get_linker_index_asymmetric( core::pose::Pose const &pose, utility::vector1< core::Size > const & res_indices ) const override;

	/// @brief Given indices of two cysteine residues that are already linked to pieces of a linker, get
	/// of the indices of the symmetric pieces of the linker.
	/// @details Throws an error if a residue is not linked to something.
	void get_linker_indices_symmetric( core::pose::Pose const &pose, utility::vector1< core::Size > const & res_indices, utility::vector1< core::Size > & linker_indices ) const override;

	/// @brief Given a pose with residues selected to be linked by a 1,4-bis(bromomethyl)benzene crosslinker,
	/// determine whether the residues are too far apart.
	/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	bool filter_by_sidechain_distance_asymmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const &filter_multiplier ) const override;

	/// @brief Given a pose with residues selected to be linked by a 1,4-bis(bromomethyl)benzene crosslinker,
	/// determine whether the residues are too far apart.  This version is for symmetric poses.
	/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	bool filter_by_sidechain_distance_symmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const &filter_multiplier ) const override;

	/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
	/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	bool filter_by_constraints_energy_asymmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const &filter_multiplier) const override;

	/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  This version is for symmetric poses.
	/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	bool filter_by_constraints_energy_symmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, bool const linker_was_added, core::Real const &filter_multiplier) const override;

	/// @brief Does this CrosslinkerMoverHelper add a residue for the linker?
	/// @details Yes, it does.
	inline bool helper_adds_linker_residue() const override { return true; }

private: // private methods

	/// @brief Confirm that the symmetry type is correct.
	/// @details Throws with an error message if is it is not.
	void confirm_symmetry( std::string const & errmsg_header ) const;

private: // data

};

} //crosslinker
} //protocols
} //cyclic_peptide

#endif //protocols/cyclic_peptide_crosslinker_1_4_BBMB_Helper_hh
