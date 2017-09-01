// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/TMA_Helper.hh
/// @brief A derived class of the CrosslinkerMoverHelper base class, used to set up
/// the trimesic acid (TMA) cross-linker.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_crosslinker_TMA_Helper_hh
#define INCLUDED_protocols_cyclic_peptide_crosslinker_TMA_Helper_hh

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/TMA_Helper.fwd.hh>
#include <protocols/cyclic_peptide/crosslinker/CrosslinkerMoverHelper.hh>

// Protocol headers

// Core headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>

// Basic/Utility headers
#include <core/types.hh>

namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

/// @brief A derived class of the CrosslinkerMoverHelper base class, used to set up
/// the trimesic acid (TMA) cross-linker.
class TMA_Helper : public CrosslinkerMoverHelper {

public: //Constructors

	/// @brief Default constructor
	TMA_Helper();

	/// @brief Copy constructor
	TMA_Helper( TMA_Helper const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~TMA_Helper() override;


public: // public methods

	/// @brief Given a pose and a selection of exactly three residues, add the TMA linker,
	/// align it crudely to the selected residues, and set up covalent bonds.
	void add_linker_asymmetric(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const override;

	/// @brief Given a pose and a linker, add bonds between the linker and the residues that coordinate the linker.
	/// @details Called by add_linker_asymmetric().  Version for asymmetric poses.
	void add_linker_bonds_asymmetric(core::pose::Pose &pose, utility::vector1< core::Size > const & res_indices, core::Size const linker_index ) const override;

	/// @brief Given a pose and a selection of exactly three residues, add the TMA linker,
	/// align it crudely to the selected residues, and set up covalent bonds.
	/// @details Version for symmetric poses.
	void add_linker_symmetric(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const override;

	/// @brief Given a pose and a TMA linker, add bonds between the TMA and the residues that coordinate the linker.
	/// @details Called by add_linker_symmetric().  Version for symmetric poses.
	void add_linker_bonds_symmetric(core::pose::Pose &pose, core::Size const res1, core::Size const linker_index1, core::Size const linker_index2 ) const override;

	/// @brief Given a selection of exactly three residues that have already been connected to a trimesic acid crosslinker,
	/// add constraints for the crosslinker.
	void add_linker_constraints_asymmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const override;

	/// @brief Given a selection of exactly three residues that have already been connected to a trimesic acid crosslinker,
	/// add constraints for the crosslinker.  This version is for symmetric poses.
	void add_linker_constraints_symmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection, bool const linker_was_added) const override;

	/// @brief Given indices of three residues that are already linked to a TMA, get the index
	/// of the TMA residue.
	/// @details Throws an error if the three residues are not all linked to the same TMA residue.
	core::Size get_linker_index_asymmetric( core::pose::Pose const &pose, utility::vector1 < core::Size > const & res_indices ) const override;

	/// @brief Given indices of three residues that are already linked to pieces of a linker, get
	/// of the indices of the symmetric pieces of the linker.
	/// @details Throws an error if a residue is not linked to something.
	void get_linker_indices_symmetric( core::pose::Pose const &pose, utility::vector1< core::Size > const & res_indices, utility::vector1< core::Size > & linker_indices ) const override;

	/// @brief Given a pose with residues selected to be linked by a trimesic acid crosslinker,
	/// determine whether the residues are too far apart.
	/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	bool filter_by_sidechain_distance_asymmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const &filter_multiplier ) const override;

	/// @brief Given a pose with residues selected to be linked by a trimesic acid crosslinker,
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

	/// @brief Optional steps that the helper can apply before every relaxation round.
	/// @details overrides default (doing nothing) to update positions of amide bond dependent atoms.
	void pre_relax_round_update_steps(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const &selection, bool const whole_structure, bool const symmetric, bool const linker_was_added) const override;

	/// @brief Optional steps that the helper can apply after every relaxation round.
	/// @details overrides default (doing nothing) to update positions of amide bond dependent atoms.
	void post_relax_round_update_steps(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const &selection, bool const whole_structure, bool const symmetric, bool const linker_was_added) const override;

private: // private methods

	/// @brief Given a residue type, return true for D- or L-Lys, Orn, DAB, or DAP/DPP, and false otherwise.
	/// @details Returns true only for the residue types that I can connect to TMA.
	bool is_allowed_residue_type( core::chemical::ResidueType const &type ) const;

	/// @brief Returns a string listing the allowed residue types that can connect to TMA.
	/// @details Intended to be human-readable.  Used for error messages.
	std::string list_allowed_residue_types() const;

	/// @brief Given a pose and a position that is one of the allowed types, add the sidechain-conjugation
	/// type or patch to the position.
	void add_sidechain_conjugation( core::pose::Pose &pose, core::Size const position ) const;

	/// @brief Given a residue type, get the name of the sidechain amide atom.
	/// @details Only works for the types allowed to connect to TMA (LYS, ORN, DAP/DPP, DAB and their D-versions).
	std::string get_sidechain_amide_name( core::chemical::ResidueType const &type ) const;

	/// @brief Given a residue type, get the name of the last carbon in the side chain.
	/// @details Only works for the types allowed to connect to TMA (LYS, ORN, DAP/DPP, DAB and their D-versions).
	std::string get_last_sidechain_carbon_name( core::chemical::ResidueType const &type ) const;

	/// @brief Given a residue in a pose, get the index of its first sidechain connection.
	///
	core::Size get_sidechain_connect_index( core::conformation::Residue const &res ) const;

	/// @brief Given one of the residue types that TMA can connect to, return the max distance between the CA
	/// and the NZ.  (Slightly padded).
	core::Real get_max_sidechain_length( core::chemical::ResidueType const &restype ) const;

	/// @brief Given a pose and three connector residue indices, create a new pose with TMA placed somewhere in space that
	/// crudely aligns it to the connector residues.
	void place_tma_asymmetric( core::pose::Pose &tma_pose, core::pose::Pose const &pose, core::Size const res1, core::Size const res2, core::Size const res3 ) const;

	/// @brief Given a pose containing only an asymmetric TMA molecule, align a symmetric TMA to
	/// the first third of the asymmetric TMA, place that in symm_pose.
	void place_tma_symmetric( core::pose::Pose &symm_pose, core::pose::Pose const &asymm_pose ) const;

	/// @brief Given a pose and a selection of three residues that connect to a TMA residue,
	/// update the positions of carbonyl oxygen and amide hydrogen atoms in the amide bonds
	/// connecting the TMA to the sidechains.
	void update_tma_amide_bond_dependent_atoms(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const &selection, bool const symmetric, bool const tma_was_added) const;

	/// @brief Given a pose and a selection of three residues that connect to a TMA residue,
	/// update the positions of carbonyl oxygen and amide hydrogen atoms in the amide bonds
	/// connecting the TMA to the sidechains.
	/// @details Version for asymmetric poses.
	void update_tma_amide_bond_dependent_atoms_asymmetric(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const &selection) const;

	/// @brief Given a pose and a selection of three residues that connect to a TMA residue,
	/// update the positions of carbonyl oxygen and amide hydrogen atoms in the amide bonds
	/// connecting the TMA to the sidechains.
	/// @details Version for symmetric poses.
	void update_tma_amide_bond_dependent_atoms_symmetric(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const &selection, bool const tma_was_added) const;


private: // data

};

} //crosslinker
} //protocols
} //cyclic_peptide

#endif //protocols/cyclic_peptide_crosslinker_TMA_Helper_hh
