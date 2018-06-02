// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/TrigonalPlanarMetal_Helper.hh
/// @brief A helper class for setting up trigonal planarly-coordinated metals
/// like Zn or Cu.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_crosslinker_TrigonalPlanarMetal_Helper_hh
#define INCLUDED_protocols_cyclic_peptide_crosslinker_TrigonalPlanarMetal_Helper_hh

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/TrigonalPlanarMetal_Helper.fwd.hh>
#include <protocols/cyclic_peptide/crosslinker/CrosslinkerMoverHelper.hh>
#include <protocols/cyclic_peptide/crosslinker/Metal_HelperBase.hh>

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

/// @brief A helper class for setting up trigonal planarly-coordinated metals
/// like Zn or Cu.
class TrigonalPlanarMetal_Helper : public Metal_HelperBase {

public: //Constructors

	/// @brief Default constructor
	TrigonalPlanarMetal_Helper( std::string const &metal_name_in = "Zn" );

	/// @brief Copy constructor
	TrigonalPlanarMetal_Helper( TrigonalPlanarMetal_Helper const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	virtual ~TrigonalPlanarMetal_Helper();


public: // public methods


protected: // protected methods


protected: // private methods

	/// @brief Check that the correct number of residues have been selected, that they are within the pose, and that they are allowed residue types.
	void check_residue_indices_valid( utility::vector1< core::Size > const &indices, core::pose::Pose const &pose ) const override;

	/// @brief Given a residue type, check whether it's an allowed residue type for trigonal planarly coordinating metals.
	/// @details Returns "true" for pass (allowed type) and "false" for failure (prohibited type).  Currently, allowed types are
	/// L- and D-histidine, L- or D-aspartate, L- or D-glutamate, the beta-3-amino acid equivalents.
	bool is_allowed_type( core::chemical::ResidueType const &type ) const override;

	/// @brief Given a pose, a list of residues, and indices i and j in that list, add angle constraints between the two residues specified.
	void add_angle_constraints( core::pose::Pose &pose, utility::vector1< core::Size > const &res_indices, core::Size const i, core::Size const j) const override;

	/// @brief Add improper dihedral constraints to keep the metal in the plane of the liganding atoms.
	void add_dihedral_constraints( core::pose::Pose &pose, utility::vector1< core::Size > const &res_indices ) const override;

	/// @brief Given a metal type and a metal-liganding atom type, return the ideal bond length.
	/// @details This method must be updated if either enum is expanded.
	core::Real const & ideal_bond_length( Metal_HelperBase_Metal const metal_type, Metal_HelperBase_MetalLigand const ligand_type ) const override;

	/// @brief Check that the symmetry type is one of a few compatible types.
	/// @details Allowed types are C3.
	void check_compatible_symmetry_type() const override;

};

} //crosslinker
} //protocols
} //cyclic_peptide

#endif //protocols/cyclic_peptide_crosslinker_TrigonalPlanarMetal_Helper_hh
