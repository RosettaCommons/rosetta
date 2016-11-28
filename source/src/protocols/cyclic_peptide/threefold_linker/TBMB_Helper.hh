// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/threefold_linker/TBMB_Helper.hh
/// @brief A derived class of the ThreefoldLinkerMoverHelper base class, used to set up
/// the 1,3,5-tris(bromomethyl)benzene (TBMB) cross-linker.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_threefold_linker_TBMB_Helper_hh
#define INCLUDED_protocols_cyclic_peptide_threefold_linker_TBMB_Helper_hh

// Unit headers
#include <protocols/cyclic_peptide/threefold_linker/TBMB_Helper.fwd.hh>
#include <protocols/cyclic_peptide/threefold_linker/ThreefoldLinkerMoverHelper.hh>

// Protocol headers

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <core/types.hh>

namespace protocols {
namespace cyclic_peptide {
namespace threefold_linker {

/// @brief A derived class of the ThreefoldLinkerMoverHelper base class, used to set up
/// the 1,3,5-tris(bromomethyl)benzene (TBMB) cross-linker.
class TBMB_Helper : public ThreefoldLinkerMoverHelper {

public: //Constructors

	/// @brief Default constructor
	TBMB_Helper();

	/// @brief Copy constructor
	TBMB_Helper( TBMB_Helper const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	virtual ~TBMB_Helper();


public: // public methods

	/// @brief Given a pose and a selection of exactly three residues, add the TBMB linker,
	/// align it crudely to the selected residues, and set up covalent bonds.
	virtual void add_linker_asymmetric(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const;

	/// @brief Given a pose and a linker, add bonds between the linker and the residues that coordinate the linker.
	/// @details Called by add_linker_asymmetric().  Version for asymmetric poses.
	virtual void add_linker_bonds_asymmetric(core::pose::Pose &pose, core::Size const res1, core::Size const res2, core::Size const res3, core::Size const linker_index ) const;

	/// @brief Given a pose and a selection of exactly three residues, add the TBMB linker,
	/// align it crudely to the selected residues, and set up covalent bonds.
	/// @details Version for symmetric poses.
	virtual void add_linker_symmetric(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const;

	/// @brief Given a selection of exactly three residues that have already been connected to a 1,3,5-tris(bromomethyl)benzene crosslinker,
	/// add constraints for the crosslinker.
	virtual void add_linker_constraints_asymmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const;

	/// @brief Given a selection of exactly three residues that have already been connected to a 1,3,5-tris(bromomethyl)benzene crosslinker,
	/// add constraints for the crosslinker.  This version is for symmetric poses.
	virtual void add_linker_constraints_symmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection, bool const linker_was_added) const;

	/// @brief Given indices of three cysteine residues that are already linked to a TBMB, get the index
	/// of the TBMB residue.
	/// @details Throws an error if the three cysteines are not all linked to the same TBMB residue.
	virtual core::Size get_linker_index_asymmetric( core::pose::Pose const &pose, core::Size const res1, core::Size const res2, core::Size const res3 ) const;

	/// @brief Given indices of three cysteine residues that are already linked to pieces of a linker, get
	/// of the indices of the symmetric pieces of the linker.
	/// @details Throws an error if a residue is not linked to something.
	virtual void get_linker_indices_symmetric( core::pose::Pose const &pose, core::Size const res1, core::Size const res2, core::Size const res3, core::Size &linker_index1, core::Size &linker_index2, core::Size &linker_index3 ) const;

	/// @brief Given a pose with residues selected to be linked by a 1,3,5-tris(bromomethyl)benzene crosslinker,
	/// determine whether the residues are too far apart.
	/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	virtual bool filter_by_sidechain_distance_asymmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const &filter_multiplier ) const;

	/// @brief Given a pose with residues selected to be linked by a 1,3,5-tris(bromomethyl)benzene crosslinker,
	/// determine whether the residues are too far apart.  This version is for symmetric poses.
	/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	virtual bool filter_by_sidechain_distance_symmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const &filter_multiplier ) const;

	/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
	/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	virtual bool filter_by_constraints_energy_asymmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const &filter_multiplier) const;

	/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  This version is for symmetric poses.
	/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	virtual bool filter_by_constraints_energy_symmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, bool const linker_was_added, core::Real const &filter_multiplier) const;

private: // private methods

	/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  Both the symmetric and asymmetric versions call this code.
	/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
	bool filter_by_constraints_energy( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, bool const symmetric, bool const linker_was_added, core::Real const &filter_multiplier ) const;

private: // data

};

} //threefold_linker
} //protocols
} //cyclic_peptide

#endif //protocols/cyclic_peptide_threefold_linker_TBMB_Helper_hh
