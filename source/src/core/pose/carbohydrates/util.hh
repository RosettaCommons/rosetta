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
/// @author  Labonte <JWLabonte@jhu.edu> and Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_core_pose_carbohydrates_util_HH
#define INCLUDED_core_pose_carbohydrates_util_HH

// Unit header
#include <core/pose/Pose.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/chemical/carbohydrates/carbohydrate_data_structures.hh>
#include <core/chemical/carbohydrates/LinkageType.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <string>


namespace core {
namespace pose {
namespace carbohydrates {

// Helper Functions
/// @brief  Use a saccharide residue's connections to find the residue from which it follows or branches.
core::uint
find_seqpos_of_saccharides_parent_residue( conformation::Residue const & residue );

/// @brief  Return pointers to the two residues of the glycosidic bond.
std::pair< conformation::ResidueCOP, conformation::ResidueCOP >
get_glycosidic_bond_residues( Pose const & pose,
	uint const sequence_position );


/// @brief  Use a saccharide residue's connections to find its linkage number on the previous residue.
core::uint
get_linkage_position_of_saccharide_residue( Pose const & pose, uint const seqpos );

core::uint
get_linkage_position_of_saccharide_residue( conformation::Residue const & rsd, conformation::Residue const & parent_rsd);

///@brief Get whether the glycosidic linkage between the residue and previous residue (parent residue) has an exocyclic carbon.
/// Does not currently work for aa->glycan.  Returns false if previous residue is not carbohydrate.
bool
has_exocyclic_glycosidic_linkage( Pose const & pose, uint const seqpos);

///@brief Get whether the glycosidic linkage between the residue and previous residue (parent residue) has an exocyclic carbon.
/// Does not currently work for aa->glycan.  Returns false if previous residue is not carbohydrate.
bool
has_exocyclic_glycosidic_linkage( conformation::Residue const & rsd, conformation::Residue const & parent_rsd );

/// @brief  Return the AtomIDs of the four phi torsion reference atoms.
utility::vector1< id::AtomID > get_reference_atoms_for_phi( Pose const & pose, uint const sequence_position );

/// @brief  Return the AtomIDs of the four psi torsion reference atoms.
utility::vector1< id::AtomID > get_reference_atoms_for_psi( Pose const & pose, uint const sequence_position );

/// @brief  Return the AtomIDs of the four omega torsion reference atoms.
utility::vector1< id::AtomID > get_reference_atoms_for_1st_omega( Pose const & pose, uint const sequence_position );

/// @brief  Return the AtomIDs of the four omega2 torsion reference atoms.
utility::vector1< id::AtomID > get_reference_atoms_for_2nd_omega( Pose const & pose, uint const sequence_position );

/// @brief  Return the AtomIDs of the four reference atoms for the requested torsion.
utility::vector1< id::AtomID > get_reference_atoms( uint const torsion_id,
	Pose const & pose,
	uint const sequence_position );


// Virtual Atom Alignment
/// @brief  Set coordinates of virtual atoms (used as angle reference points) within a saccharide residue of the given
/// conformation.
void align_virtual_atoms_in_carbohydrate_residue( conformation::Conformation & conf, uint const sequence_position );

/// @brief  Set coordinates of virtual atoms (used as angle reference points) within a saccharide residue of the given
/// pose.
void align_virtual_atoms_in_carbohydrate_residue( Pose & pose, uint const sequence_position );


// TorsionID Queries
/// @brief  Is this is the phi torsion angle of a glycosidic linkage?
bool is_glycosidic_phi_torsion( Pose const & pose, id::TorsionID const & torsion_id );

/// @brief  Is this is the psi torsion angle of a glycosidic linkage?
bool is_glycosidic_psi_torsion( Pose const & pose, id::TorsionID const & torsion_id );

/// @brief  Is this is an omega torsion angle of a glycosidic linkage?
bool is_glycosidic_omega_torsion( Pose const & pose, id::TorsionID const & torsion_id );


// Torsion Access
// Getters

///@brief Get the number of glycosidic torsions for this residue.  Up to 4 (omega2).
Size get_n_glycosidic_torsions_in_res(
	Pose & pose,
	uint const sequence_position);
	
/// @brief  Return the requested torsion angle between a saccharide residue of the given pose and the previous residue.
core::Angle get_glycosidic_torsion( uint const torsion_id, Pose const & pose, uint const sequence_position );


// Setters
/// @brief  Set the requested torsion angle between a saccharide residue of the given pose and the previous residue.
void set_glycosidic_torsion(
	uint const torsion_id,
	Pose & pose,
	uint const sequence_position,
	core::Angle const setting );

// Glycosylation
/// @brief  Idealize the glycosidic torsion angles for the last n glycan residues added or built.
void idealize_last_n_glycans_in_pose( Pose & pose, Size const n_glycans_added );

/// @brief  Glycosylate the Pose at the given sequence position and atom using an IUPAC sequence.
void glycosylate_pose(
	Pose & pose,
	uint const sequence_position,
	std::string const & atom_name,
	std::string const & iupac_sequence,
	bool const idealize_linkages = true );

/// @brief  Glycosylate the Pose at the given sequence position using an IUPAC sequence.
void glycosylate_pose(
	Pose & pose,
	uint const sequence_position,
	std::string const & iupac_sequence,
	bool const idealize_linkages = true );

/// @brief  Glycosylate the Pose at the given sequence position and atom using a .GWS or IUPAC sequence file.
void glycosylate_pose_by_file(
	Pose & pose,
	uint const sequence_position,
	std::string const & atom_name,
	std::string const & filename,
	bool const idealize_linkages = true );

/// @brief  Glycosylate the Pose at the given sequence position using a .GWS or IUPAC sequence file.
void glycosylate_pose_by_file(
	Pose & pose,
	uint const sequence_position,
	std::string const & filename,
	bool const idealize_linkages = true );



/// @brief  Set the dihedral angles involved in a glycosidic linkage based on statistical data.
void
set_dihedrals_from_linkage_conformer_data( Pose & pose,
	uint const upper_residue,
	core::chemical::carbohydrates::LinkageConformerData const & conformer,
	bool idealize = true,
	bool use_prob_for_sd = false );


/// @brief Get the linkage type for a particular residue.
/// Be warned: If NO STATISTICS for that residue exist for SugarBB/CHI (such as not a pyranose, will return LINKAGE_NA )
chemical::carbohydrates::LinkageType
get_linkage_type_for_residue_for_CHI( core::Size torsion_id, conformation::Residue const & rsd,
	pose::Pose const & pose);

utility::vector1< chemical::carbohydrates::LinkageType >
get_linkage_types_for_dihedral( core::Size torsion_id );


/// @brief Remove ALL/Any branch points from a carbohydrate or aa residue.
void
remove_carbohydrate_branch_point_variants( Pose & pose, core::Size const seqpos );





/////////////////////////////////////////// Branch Deletion ////////////////////////////////////////////
///
///



/// @brief Delete the glycan from this residue onward toward the end of the branch.  Like chopping off a tree trunk at position resnum (not including the resnum). Also known as defoliating.
///  If resnum is the protein branch point, will change variant.
//   If no more carbohydrates exist in the pose, will change the pose status.
void
delete_carbohydrate_branch(
	Pose & pose,
	uint const delete_to);


/// @brief Get residues further down the branch from this residue.  starting_position ->
///  May require a better name.
///  Returns pair of all_upstream_residues, tips.
///  Tips are the ends of linear glycan branches.

std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > >
get_carbohydrate_residues_upstream(
	Pose const & pose,
	uint const starting_position);


/// @brief  Recursive function to get branches of a set of residues, etc.
///  list_of_residues and tips are arrays are non-const references and modified by this function.
///
///  Children Residues:  Residue nums of parent residue connected that we are interested in finding connected branchs.
///  List Of  Residues:  All the residue nums of the branching from children residues
///  Tips:  All 'ends' of all the branches found using this function.
///
///  See Also: get_carbohydrate_residues_upstream
///            trim_carbohydrate_branch_from_X
///
void
get_branching_residues( Pose const & pose,
	Size parent_residue,
	utility::vector1< Size > & children_residues,
	utility::vector1< Size > & list_of_residues,
	utility::vector1< Size > & tips );

/// @brief  Find all children residues, list of residues, and any found tips from a given residue not including parent
///
///  Children Residues:  Filled in list of children residues found if not tips.
///  List Of  Residues:  All the residue nums found.
///  Tips:  All 'ends' of of children found.
///
///  See Also: get_carbohydrate_residues_upstream
///            trim_carbohydrate_branch_from_X
///
void
fill_upstream_children_res_and_tips( Pose const & pose,
	Size res,
	Size parent_residue,
	utility::vector1< Size > & children_residues,
	utility::vector1< Size > & list_of_residues,
	utility::vector1< Size > & tips );


/// @brief Get all residue numbers in order from the tip to (and not including) stop_at_residue or a branch point.
///  All residue numbers are the tip or a linear polymer of glycans.
///  Useful for glycan stripping.
///
utility::vector1< core::Size >
get_resnums_in_leaf( Pose const & pose, Size tip_residue, Size stop_at_residue);


/// @brief Delete a leaf of glycan residues.  Use the ReferencePose to associate residue numbers.
/// @details
///  Uses delete_residue_slow as the debug_assert in delete_polymer_residue causes a crash.
///  This is a bug in the FoldTree where chemical edges are being treated as Jumps.
///  This is being addressed in a branch and a FoldTree.
void
delete_leaf( Pose & pose, utility::vector1< Size > leaf, std::string ref_pose_name = "temp_ref_pose" );


}  // namespace carbohydrates
}  // namespace pose
}  // namespace core

#endif  // INCLUDED_core_pose_carbohydrates_util_HH
