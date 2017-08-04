// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/util/metalloproteins_util.hh
/// @brief  Utilities for working with metalloproteins
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_util_metalloproteins_util_hh
#define INCLUDED_core_util_metalloproteins_util_hh

// Project headers
#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>

// Package headers
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>
#include <core/pose/Pose.fwd.hh>
#include <map>
#include <set>
namespace core {
namespace util {

#define DEFAULT_DIST_CUTOFF_MULTIPLIER 1.05 //By default, the distance cutoff multiplier gives a cutoff radius just slightly larger than the sum of the atoms' Lennard-Jones radii

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief: Adds an arbitrary covalent linkage between two atoms (resA_At and resB_At) in two residues (at positions resA_pos and resB_pos).
/// @details:  This is useful for adding covalent linkages between metal-binding side-chains and metal atoms.  This code was shamelessly
/// stolen from Florian's EnzConstraintParameters.cc in protocols/toolbox/match_enzdes_utils, and was modified to permit deletion of
/// unnecessary protons.  NOTE: THIS CODE MODIFIES THE RESIDUE TYPE LIST, AND IS CURRENTLY NOT THREADSAFE.
/// @author:  Vikram K. Mulligan (vmullig@uw.edu), Florian Richter (flosopher@gmail.com)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
add_covalent_linkage(
	core::pose::Pose & pose,
	Size const resA_pos,
	Size const resB_pos,
	Size const resA_At,
	Size const resB_At,
	bool const remove_hydrogens //Should extraneous hydrogens on the bonding atoms be removed?
);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///@brief: Set of element strings used to identify metal ions.
std::set< std::string > const METALS = { "NA", "K", "MG", "LI", "CA", "CU", "CO", "NI", "FE", "ZN", "CR", "MN", "MO", "PD", "PT", "AG", "HG", "LA", "FE2"};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///@brief Function that can find contacts to metal ions that are part of a larger complex.
///@details Identifies metal ions based on the METALS set defined in this file and passes
///those atoms to find_metalbinding_atoms_helper.
///@author Sharon Guffy (guffy@email.unc.edu )
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::map< core::Size, utility::vector1< core::id::AtomID > >
find_metalbinding_atoms_for_complex(
	core::pose::Pose const &pose,
	core::Size const metal_position,
	core::Real const dist_cutoff_multiplier
);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///@brief Helper function that finds contacts for metal atoms specified by AtomID
///@details Code modified from Vikram's code originally in find_metalbinding_atoms
///@author Sharon Guffy (guffy@email.unc.edu )
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< core::id::AtomID >
find_metalbinding_atoms_helper(
	core::pose::Pose const &pose,
	core::id::AtomID const & metal_atom,
	core::Real const dist_cutoff_multiplier
);






/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function that generates a list of metal-binding atoms that coordinate a metal in a protein.
/// @details This function generates the list by looping through all residues and checking all metal-binding atoms of all
/// metal-binding residues, so it's not super speedy.
/// Inputs:
///  pose (The pose that we'll operate on, unchanged by operation)
///  metal_postion (The residue number of the metal)
///  dist_cutoff_multiplier (A float for the distance cutoff multiplier; the cutoff is the sum of the Lennard-Jones radii times the multiplier)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1 < core::id::AtomID >
find_metalbinding_atoms(
	core::pose::Pose const &pose,
	core::Size const metal_position,
	core::Real const dist_cutoff_multiplier = DEFAULT_DIST_CUTOFF_MULTIPLIER
);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to add covalent linkages between a metal atom and all the liganding atoms provided in a vector of AtomIDs.
/// @details The inputs are:
///    pose (The pose to be modified)
///    metal_position (The residue number of the metal in the pose)
///    liganding_atomids (A list of AtomIDs on other residues that will be covalently linked to the metal)
///    remove_hydrogens (Should hydrogens on the liganding atoms be removed automatically?  Default true.)
/// This function uses core::pose::add_covalent_linkage, which can strip off extraneous hydrogens.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
add_covalent_linkages_to_metal(
	core::pose::Pose &pose,
	core::Size const metal_position,
	utility::vector1 < core::id::AtomID > &liganding_atomids,
	bool const remove_hydrogens = true
);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to auto-detect and add covalent connections to all metal ions in a pose.
/// @details This function iteratively calls auto_setup_metal_bonds.
/// Inputs:
///  pose (The pose that we'll operate on, changed by operation)
///  dist_cutoff_multiplier (A float for the distance cutoff multiplier; the cutoff is the sum of the Lennard-Jones radii times the multiplier)
///   remove_hydrogesn (Should hydrogens on the liganding atoms be auto-removed?  Default true.)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
auto_setup_all_metal_bonds (
	core::pose::Pose &pose,
	core::Real const dist_cutoff_multiplier = DEFAULT_DIST_CUTOFF_MULTIPLIER,
	bool const remove_hydrogens = true
);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to set up distance and angle constraints for a specified metal ion.
/// @details Modified/taken from auto_setup_all_metal_constraints to allow constraint setup for only specified metal residues.
/// Does not modify score function weights.
/// Inputs:
///  pose (The pose that we'll operate on, changed by operation)
///  metal_position (The residue number of the metal to be constrained)
///  distance_constraint_multiplier (A float for the strength of the metal - binding atom distance constraint.  A value of 2.0 doubles
///   it, for example.)
///  angle_constraint_multiplier (A float for the strength of the metal - binding atom - binding atom parent angle constraint.)
/// @author Sharon Guffy (guffy@email.unc.edu); originally written by Vikram K. Mulligan (vmulligan@uw.edu) as part of auto_setup_all_metal_constraints
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
add_constraints_to_metal(
	core::pose::Pose &pose,
	core::Size const metal_position,
	core::Real const distance_constraint_multiplier,
	core::Real const angle_constraint_multiplier
);


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to set up distance and angle constraints between metals and the residues that bind them.
/// @details This function constrains the distances to be whatever they are in the input pose.  This version
/// does not set the weights for the constraints terms in the scorefunction.
/// Inputs:
///  pose (The pose that we'll operate on, changed by operation)
///   distance_constraint_multiplier (A float for the strength of the metal - binding atom distance constraint.  A value of 2.0 doubles
///   it, for example.)
///   angle_constraint_multiplier (A float for the strength of the metal - binding atom - binding atom parent angle constraint.)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
auto_setup_all_metal_constraints(
	core::pose::Pose &pose,
	core::Real const distance_constraint_multiplier,
	core::Real const angle_constraint_multiplier
);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Function to set up distance and angle constraints between metals and the residues that bind them.
/// @details This function constrains the distances to be whatever they are in the input pose.  This version
/// sets the weights for the constraints terms in the scorefunction to 1.0 if they're off, or scales the
/// constraints themselves appropriately if they're already on.
/// Inputs:
///  pose (The pose that we'll operate on, changed by operation)
///   sfxn (An owning pointer to the scorefunction, changed by operation)
///   distance_constraint_multiplier (A float for the strength of the metal - binding atom distance constraint.  A value of 2.0 doubles
///   it, for example.)
///   angle_constraint_multiplier (A float for the strength of the metal - binding atom - binding atom parent angle constraint.)
/// @author Vikram K. Mulligan (vmulligan@uw.edu)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
auto_setup_all_metal_constraints(
	core::pose::Pose &pose,
	core::scoring::ScoreFunctionOP sfxn,
	core::Real const distance_constraint_multiplier,
	core::Real const angle_constraint_multiplie
);

} // util
} // core

#endif // INCLUDED_core_util_metalloproteins_util_hh
