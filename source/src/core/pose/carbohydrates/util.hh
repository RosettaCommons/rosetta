// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/pose/carbohydrates/util.hh
/// @brief   Utility function declarations for carbohydrate-containing poses.
/// @author  labonte

#ifndef INCLUDED_core_pose_carbohydrates_util_HH
#define INCLUDED_core_pose_carbohydrates_util_HH

// Unit header
#include <core/pose/Pose.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <string>


namespace core {
namespace pose {
namespace carbohydrates {

// Helper Functions
/// @brief  Scan through a saccharide residue's connections to find the residue from which it follows or branches.
core::uint find_seqpos_of_parent_residue(conformation::Residue const & residue);

/// @brief  Return pointers to the two residues of the glycosidic bond.
std::pair<conformation::ResidueCOP, conformation::ResidueCOP> get_glycosidic_bond_residues(Pose const & pose,
		uint const sequence_position);

/// @brief  Scan through the list of atoms connected to a given "connect_atom" and return the first heavy atom found.
core::uint atom_next_to_connect_atom(conformation::Residue const & residue, core::uint const connect_atom_index);


/// @brief  Return the AtomIDs of the four phi torsion reference atoms.
utility::vector1<id::AtomID> get_reference_atoms_for_phi(Pose const & pose, uint const sequence_position);

/// @brief  Return the AtomIDs of the four psi torsion reference atoms.
utility::vector1<id::AtomID> get_reference_atoms_for_psi(Pose const & pose, uint const sequence_position);

/// @brief  Return the AtomIDs of the four omega torsion reference atoms.
utility::vector1<id::AtomID> get_reference_atoms_for_1st_omega(Pose const & pose, uint const sequence_position);

/// @brief  Return the AtomIDs of the four reference atoms for the requested torsion.
utility::vector1<id::AtomID> get_reference_atoms(uint const torsion_id,
		Pose const & pose,
		uint const sequence_position);


// Virtual Atom Alignment
/// @brief  Set coordinates of virtual atoms (used as angle reference points) within a saccharide residue of the given
/// conformation.
void align_virtual_atoms_in_carbohydrate_residue(conformation::Conformation & conf, uint const sequence_position);

/// @brief  Set coordinates of virtual atoms (used as angle reference points) within a saccharide residue of the given
/// pose.
void align_virtual_atoms_in_carbohydrate_residue(Pose & pose, uint const sequence_position);


// Torsion Access
// Getters
/// @brief  Return the requested torsion angle between a saccharide residue of the given pose and the previous residue.
core::Angle get_glycosidic_torsion(uint const torsion_id, Pose const & pose, uint const sequence_position);


// Setters
/// @brief  Set the requested torsion angle between a saccharide residue of the given pose and the previous residue.
void set_glycosidic_torsion(uint const torsion_id,
		Pose & pose,
		uint const sequence_position,
		core::Angle const setting);

}  // namespace carbohydrates
}  // namespace pose
}  // namespace core

#endif  // INCLUDED_core_pose_carbohydrates_util_HH
