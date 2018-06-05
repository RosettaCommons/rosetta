// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/Metal_HelperBase.hh
/// @brief Headers for a base class for setting up metals.  This is a pure virtual class that must be subclassed for
/// specific metal geometries.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_crosslinker_Metal_HelperBase_hh
#define INCLUDED_protocols_cyclic_peptide_crosslinker_Metal_HelperBase_hh

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/Metal_HelperBase.fwd.hh>
#include <protocols/cyclic_peptide/crosslinker/CrosslinkerMoverHelper.hh>

// Protocol headers

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

/// @brief Metal types compatible with this helper.
/// @details If you add to this list, you need to add entries to
/// the switch statement in metal_type_string_from_enum(), as well as
/// the map in ideal_bond_length().
enum Metal_HelperBase_Metal {
	MH_Zn = 1, //Keep this first.
	MH_Fe2,
	MH_Ni2,
	MH_unknown_metal, //Keep this second-to-last.
	MH_end_of_list = MH_unknown_metal //Keep this last.
};

/// @brief Metal-liganding atom types compatible with this helper.
/// @details If you add to this list, you need to add entries to
/// the map in ideal_bond_length().
enum Metal_HelperBase_MetalLigand {
	MHLigand_Nd_histidine = 1, //Keep first.
	MHLigand_Ne_histidine,
	MHLigand_O_carboxyl,
	MHLigand_N_pyridine,
	MHLigand_S_cysteine,
	MHLigand_unknown_ligand, //Keep second-to-last.
	MHLigand_end_of_list = MHLigand_unknown_ligand //Keep last.
};

/// @brief A base class for setting up metals.  This is a pure virtual class that must be subclassed for
/// specific metal geometries..
class Metal_HelperBase : public CrosslinkerMoverHelper {

public: //Constructors

	/// @brief Default constructor
	Metal_HelperBase( std::string const & metal_name_in = "Zn" );

	/// @brief Copy constructor
	Metal_HelperBase( Metal_HelperBase const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	virtual ~Metal_HelperBase();


public: // public methods

	/// @brief Given a pose and a selection of exactly three residues, add the linker,
	/// align it crudely to the selected residues, and set up covalent bonds.
	/// @details Must be defined by derived classes.  Version for asymmetric poses.
	void add_linker_asymmetric(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const override;

	/// @brief Given a pose and a linker, add bonds between the linker and the residues that coordinate the linker.
	/// @details Can be called by add_linker_asymmetric().  Must be defined by derived classes (pure virtual).  Version for asymmetric poses.
	void add_linker_bonds_asymmetric(core::pose::Pose &pose, utility::vector1< core::Size > const & res_indices, core::Size const linker_index ) const override;

	/// @brief Given a pose and a selection of exactly three residues, add the linker,
	/// align it crudely to the selected residues, and set up covalent bonds.
	/// @details Must be defined by derived classes.  Version for symmetric poses.
	void add_linker_symmetric(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const override;

	/// @brief Given a pose and a linker, add bonds between the linker and the residues that coordinate the linker.
	/// @details Can be called by add_linker_symmetric().  Must be defined by derived classes (pure virtual).  Version for symmetric poses.
	void add_linker_bonds_symmetric(core::pose::Pose &pose, core::Size const res1, core::Size const linker_index1, core::Size const linker_index2 ) const override;

	/// @brief Given a selection of exactly three residues that have already been connected to a crosslinker,
	/// add constraints for the crosslinker.
	/// @details Must be defined by derived classes.  Version for asymmetric poses.
	void add_linker_constraints_asymmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const override;

	/// @brief Given a selection of exactly three residues that have already been connected to a crosslinker,
	/// add constraints for the crosslinker.
	/// @details Must be defined by derived classes.  Version for symmetric poses.
	void add_linker_constraints_symmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection, bool const linker_was_added) const override;

	/// @brief Given indices of four residues that are already linked to a linker, get the index
	/// of the linker.
	/// @details Not applicable for this particular crosslinker.  A "GNDN" function -- goes nowhere, does nothing.
	/// Only here because the base class has a pure virtual of this name.
	core::Size get_linker_index_asymmetric( core::pose::Pose const &pose, utility::vector1< core::Size > const & res_indices ) const override;

	/// @brief Given indices of three cysteine residues that are already linked to pieces of a linker, get
	/// of the indices of the symmetric pieces of the linker.
	/// @details Throws an error if a residue is not linked to something.  Must be defined by derived classes.
	void get_linker_indices_symmetric( core::pose::Pose const &pose, utility::vector1< core::Size > const & res_indices, utility::vector1< core::Size > & linker_indices ) const override;

	/// @brief Given a pose with residues selected to be linked by a linker, determine whether the residues are too far apart.
	/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	bool filter_by_sidechain_distance_asymmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const & filter_multiplier ) const override;

	/// @brief Given a pose with residues selected to be linked by a linker, determine whether the residues are too far apart.
	/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.  This version is for symmetric poses.
	/// @note Higher values of the filter multiplier make it more permissive.
	bool filter_by_sidechain_distance_symmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const & filter_multiplier ) const override;

	/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
	/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	bool filter_by_constraints_energy_asymmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const & filter_multiplier ) const override;

	/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  This version is for symmetric poses.
	/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	bool filter_by_constraints_energy_symmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, bool const linker_was_added, core::Real const & filter_multiplier ) const override;

	/// @brief Does this CrosslinkerMoverHelper add a residue for the linker?
	/// @details No, it does not.
	bool helper_adds_linker_residue() const override { return false; }

public: // public defined methods

	/// @brief Given a metal name, set the type.
	void set_metal_type_from_name( std::string const & name_in );

	/// @brief Get the current metal type, as a string.
	std::string metal_type_string() const;

	/// @brief Does this helper add a residue to the pose?
	/// @details True for most crosslinkers, but overridden here to return "false" because this helper acts by
	/// patching existing residues instead of appending geometry.
	bool adds_crosslinker_residue() const override { return false; }

protected: // private methods

	/// @brief Check that the correct number of residues have been selected, that they are within the pose, and that they are allowed residue types.
	virtual void check_residue_indices_valid( utility::vector1< core::Size > const &indices, core::pose::Pose const &pose ) const = 0;

	/// @brief Given a residue type, check whether it's an allowed residue type for tetrahedrally coordinating metals.
	/// @details Returns "true" for pass (allowed type) and "false" for failure (prohibited type).  Currently, allowed types are L- and D-histidine,
	/// L- or D-aspartate, L- or D-glutamate, L- or D-cysteine, L- or D-homocysteine, and the beta-3-amino acid equivalents.
	virtual bool is_allowed_type( core::chemical::ResidueType const &type ) const = 0;

	/// @brief Given a pose, a list of residues, and indices i and j in that list, add angle constraints between the two residues specified.
	virtual void add_angle_constraints( core::pose::Pose &pose, utility::vector1< core::Size > const &res_indices, core::Size const i, core::Size const j) const = 0;

	/// @brief Given a pose and a list of residues, add dihedral constraints (e.g. improper dihedrals to enforce planarity).
	/// @details Defaults to a function that does nothing.  Can be overridden by derived classes.
	virtual void add_dihedral_constraints( core::pose::Pose &pose, utility::vector1< core::Size > const &res_indices ) const;

	/// @brief Given a metal type and a metal-liganding atom type, return the ideal bond length.
	/// @details This method must be updated if either enum is expanded.
	virtual core::Real const & ideal_bond_length( Metal_HelperBase_Metal const metal_type, Metal_HelperBase_MetalLigand const ligand_type ) const = 0;

	/// @brief Given a pose and a residue with the VIRTUAL_METAL_CONJUGATION variant type already added,
	/// set the metal-metal liganding atom bond length appropriately for the metal in question.
	/// @details Calls ideal_bond_length().
	void set_metal_bond_length( core::pose::Pose &pose, core::Size const res_index ) const;

	/// @brief Given a ResidueType with the VIRTUAL_METAL_CONJUGATION variant type already added, get the metal-liganding
	/// atom enum.
	Metal_HelperBase_MetalLigand liganding_atom_from_restype( core::chemical::ResidueType const &restype ) const;

	/// @brief Given a metal type enum, return the corresponding string.
	std::string metal_type_string_from_enum( Metal_HelperBase_Metal const metal_type ) const;

	/// @brief Given a metal type string, return the corresponding enum.
	/// @details Returns MH_unknown_metal if the string is unrecognized.
	Metal_HelperBase_Metal metal_type_enum_from_string( std::string const & metal_type_string ) const;

	// @brief Get the metal type.
	//inline Metal_HelperBase_Metal metal_type() const { return metal_type_; }

	/// @brief Check that the symmetry type is one of a few compatible types.
	virtual void check_compatible_symmetry_type() const = 0;

private: // data

	/// @brief The type of metal to add (Zn, Cu, etc.).
	/// @details Determines ideal bond lengths for metal bonds.
	Metal_HelperBase_Metal metal_type_;

};

} //crosslinker
} //protocols
} //cyclic_peptide

#endif //protocols/cyclic_peptide_crosslinker_Metal_HelperBase_hh
