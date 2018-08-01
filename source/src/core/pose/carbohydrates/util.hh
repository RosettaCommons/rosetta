// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/pose/carbohydrates/util.hh
/// @brief   Utility function declarations for carbohydrate-containing poses.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @author  Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_pose_carbohydrates_util_HH
#define INCLUDED_core_pose_carbohydrates_util_HH

// Unit Header
#include <core/pose/Pose.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/chemical/carbohydrates/LinkageConformers.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/id/SequenceMapping.hh>

// Utility Header
#include <utility/vector1.hh>

// C++ Header
#include <string>
#include <utility>


namespace core {
namespace pose {
namespace carbohydrates {

////////////////////////   Pose-Required Kinematics    ///////////////////////////////////
///
///
///

///@brief
///
///  Get a pre-configured movemap from a residue selector.
///  Use the ReturnResidueSubsetSelector to obtain from a subset.
///
///  The Rosetta Movemap is VERY different from IUPAC-designated torsions for carbohydrates.
///  NEVER attempt to create a MoveMap for carbohydrates unless you know what you are doing.
///
///@details
///
/// This will create a Movemap from the residue selector for ALL residues within it.
/// including non-carbohydrates. This includes Chemical Edge Branch points, Mainchains, ASN->glycan linakge, etc.
///
/// include_iupac_chi:
///    Include the 'carbohydrate 'side-chains' (rotatable OH groups) and any selected non-carbohydrate side-chain
///
/// include_ring_torsions:
///    Include moveable ring torsions
///
/// include_bb_torsions:
///    Include BB torsions (normal for protein) or
///    between both carboydrate or ASN/Carbohydrate residues as defined by IUPAC.
///     IE for Carbohdrate linkage 1->5, the torsions of residue 5
///
kinematics::MoveMapOP
create_glycan_movemap_from_residue_selector(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSelectorCOP selector,
	bool include_iupac_chi = true,
	bool include_glycan_ring_torsions = true,
	bool include_bb_torsions = true,
	bool cartesian=false);


///@brief Turn on/off IUPAC CHIs for a particular residue number.
void
set_glycan_iupac_chi_torsions(
	core::pose::Pose const & pose,
	core::kinematics::MoveMap & movemap,
	core::Size resnum,
	bool action=true,
	bool cartesian=false
);

///@brief Turn on/off IUPAC BBs for a particular residue number.
void
set_glycan_iupac_bb_torsions(
	core::pose::Pose const & pose,
	core::kinematics::MoveMap & movemap,
	core::Size resnum,
	bool action=true,
	bool cartesian=false);



/////////////////////////////////   Glycosylation    ////////////////////////////////////////

/// @brief  Idealize the glycosidic torsion angles for the last n glycan residues added or built.
void idealize_last_n_glycans_in_pose( Pose & pose, Size const n_glycans_added );

/// @brief  Glycosylate the Pose at the given sequence position and atom using an IUPAC sequence.
void glycosylate_pose(
	Pose & pose,
	uint const sequence_position,
	std::string const & atom_name,
	std::string const & iupac_sequence,
	bool idealize_linkages = true,
	bool keep_pdbinfo = true );

/// @brief  Glycosylate the Pose at the given sequence position using an IUPAC sequence.
void glycosylate_pose(
	Pose & pose,
	uint const sequence_position,
	std::string const & iupac_sequence,
	bool const idealize_linkages = true,
	bool keep_pdbinfo = true );

/// @brief  Glycosylate the Pose at the given sequence position and atom using a .GWS or IUPAC sequence file.
void glycosylate_pose_by_file(
	Pose & pose,
	uint const sequence_position,
	std::string const & atom_name,
	std::string const & filename,
	bool const idealize_linkages = true,
	bool keep_pdbinfo = true );

/// @brief  Glycosylate the Pose at the given sequence position using a .GWS or IUPAC sequence file.
void glycosylate_pose_by_file(
	Pose & pose,
	uint const sequence_position,
	std::string const & filename,
	bool const idealize_linkages = true,
	bool keep_pdbinfo = true );


/// @brief  Set the dihedral angles involved in a glycosidic linkage based on statistical data.
void set_dihedrals_from_linkage_conformer_data( Pose & pose,
	uint const upper_residue,
	core::chemical::carbohydrates::LinkageConformerData const & conformer,
	bool idealize = true,
	bool use_prob_for_sd = false );


///////////////////////////   On-The-Fly Helper Functions    /////////////////////////////////////

/// @brief  Return pointers to the two residues of the glycosidic bond.
std::pair< conformation::ResidueCOP, conformation::ResidueCOP > get_glycosidic_bond_residues( Pose const & pose,
	uint const sequence_position );


/// @brief  Return the AtomIDs of the four phi torsion reference atoms.
utility::vector1< id::AtomID > get_reference_atoms_for_phi( Pose const & pose, uint const sequence_position );

/// @brief  Return the AtomIDs of the four psi torsion reference atoms.
utility::vector1< id::AtomID > get_reference_atoms_for_psi( Pose const & pose, uint const sequence_position );

/// @brief  Return the AtomIDs of the four omega torsion reference atoms.
utility::vector1< id::AtomID > get_reference_atoms_for_1st_omega( Pose const & pose, uint const sequence_position );

/// @brief  Return the AtomIDs of the four omega2 torsion reference atoms.
utility::vector1< id::AtomID > get_reference_atoms_for_2nd_omega( Pose const & pose, uint const sequence_position );

/// @brief  Return the AtomIDs of the four reference atoms for the requested torsion.
utility::vector1< id::AtomID > get_reference_atoms( uint const named_torsion,
	Pose const & pose,
	uint const sequence_position );


/////////////////////////    On-The-Fly Virtual Atom Alignment    /////////////////////////////

/// @brief  Set coordinates of virtual atoms (used as angle reference points) within a saccharide residue of the given
/// pose.
void align_virtual_atoms_in_carbohydrate_residue( Pose & pose, uint const sequence_position );


//////////////////////////    On-The-Fly TorsionID Queries    //////////////////////////////////

///@brief Get the torsion num that this torsionID moves.
/// Returns 0 if not any of the 4 canonical BB torsions up to omega2.
core::Size
which_glycosidic_torsion(Pose const & pose, id::TorsionID const & torsion_id );


/// @brief  Is this is the phi torsion angle of a glycosidic linkage?
bool is_glycosidic_phi_torsion( Pose const & pose, id::TorsionID const & torsion_id );

/// @brief  Is this is the psi torsion angle of a glycosidic linkage?
bool is_glycosidic_psi_torsion( Pose const & pose, id::TorsionID const & torsion_id );

/// @brief  Is this is a 1st omega torsion angle of a glycosidic linkage?
bool is_glycosidic_omega_torsion( Pose const & pose, id::TorsionID const & torsion_id );

/// @brief  Is this is a 2nd omega torsion angle of a glycosidic linkage?
bool is_glycosidic_omega2_torsion( Pose const & pose, id::TorsionID const & torsion_id );

/// @brief  Is this is a 3rd omega torsion angle of a glycosidic linkage?
bool is_glycosidic_omega3_torsion( Pose const & pose, id::TorsionID const & torsion_id );

///@brief Base function to reduce code-duplication in torsion queries.
bool is_glycosidic_torsion(
	Pose const & pose, id::TorsionID const & torsion_id,
	core::id::MainchainTorsionType const & torsion_type );


/// @brief  Return the sequence position of the immediate downstream (child) residue affected by this torsion.
core::uint get_downstream_residue_that_this_torsion_moves( Pose const & pose, id::TorsionID const & torsion_id );


///////////////////////////    On-The-Fly Torsion Access    ///////////////////////////////////

// Getters
/// @brief Get the number of glycosidic torsions for this residue.  Up to 4 (omega2).
Size
get_n_glycosidic_torsions_in_res( Pose const & pose, uint const sequence_position );

/// @brief  Return the requested torsion angle between a saccharide residue of the given pose and the previous residue.
core::Angle get_glycosidic_torsion( uint const named_torsion, Pose const & pose, uint const sequence_position );


// Setters
/// @brief  Set the requested torsion angle between a saccharide residue of the given pose and the previous residue.
void set_glycosidic_torsion(
	uint const named_torsion,
	Pose & pose,
	uint const sequence_position,
	core::Angle const setting );


////////////////////////    Branch Information    ///////////////////////////////////

/// @brief Get residues further down the branch from this residue.  starting_position ->
/// @details Convenience function. Calls get_carbohydrate_residues_and_tips_of_branch
utility::vector1< core::Size > get_carbohydrate_residues_of_branch(
	Pose const & pose,
	uint const starting_position);


/// @brief Get tips (end residue of linear components of branches) further down the branch from this residue.  starting_position ->
/// @details Convenience function. Calls get_carbohydrate_residues_and_tips_of_branch
utility::vector1< core::Size > get_carbohydrate_tips_of_branch(
	Pose const & pose,
	uint const starting_position);


/// @brief Get residues further down the branch from this residue.  starting_position ->
///  Returns pair of all_upstream_residues, tips.
///  Tips are the ends of linear glycan branches.
std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > >
get_carbohydrate_residues_and_tips_of_branch(
	Pose const & pose,
	uint const starting_position);


/// @brief Get the carbohydrate residue connecting the protein branch point.
/// Returns 0 if branch point is not connected to carbohydrate downstream.
core::Size
get_glycan_connecting_protein_branch_point(
	Pose const & pose,
	core::Size const protein_branch_point_resnum);


///@brief Get the particular resnum from a glycan position, givin the protein branch point.
/// The glycan_position is numbered 1 -> length of glycan. This is useful for easily identifying a particular glycan position.
/// Returns 0 if that glycan_position is not part of the glycan we are interested in or not in pose.
core::Size
get_resnum_from_glycan_position(
	Pose const & pose,
	core::Size const first_glycan_resnum,
	core::Size const glycan_position);


///@brief Get the particular resnum from a glycan position, givin the protein branch point.
/// The glycan_position is numbered 1 -> length of glycan. This is useful for easily identifying a particular glycan position.
/// Returns 0 if that glycan_position is not part of the glycan we are interested in or not in pose.
core::Size
get_glycan_position_from_resnum(
	Pose const & pose,
	core::Size const first_glycan_resnum,
	core::Size const resnum);


///@brief Get all resnums from specific glycan positions given the particular glycan.  Positions correspond to a number 1 -> that corresponds to each residue.
/// Use glycan_info app for more.
select::residue_selector::ResidueSubset
get_resnums_from_glycan_positions(
	Pose const & pose,
	core::Size const first_glycan_resnum,
	utility::vector1< core::Size > const & glycan_positions);

///@brief Get all resnums for all glycans in the pose.  Positions correspond to a number 1 -> that corresponds to each residue.
/// Use glycan_info app for more. Glycan positions do not need to exist for all glcyans.
select::residue_selector::ResidueSubset
get_resnums_from_glycan_positions(
	Pose const & pose,
	utility::vector1< core::Size > const & glycan_positions);


////////////////////////////   Branch Deletion    ////////////////////////////////////////


utility::vector1< bool >
get_mainchain_children( Pose const & pose, core::Size starting_resnum, bool include_starting_resnum = false);


/// @brief Delete the glycan from this residue onward toward the end of the branch.  Like chopping off a tree trunk at position resnum (not including the resnum). Also known as defoliating.
///  If resnum is the protein branch point, will change variant.
//   If no more carbohydrates exist in the pose, will change the pose status.
void
delete_carbohydrate_branch( Pose & pose, uint const delete_to );


/// @brief  Recursive function to get branches of a set of residues, etc.
///  list_of_residues and tips are arrays are non-const references and modified by this function.
///
///  Children Residues:  Residue nums of parent residue connected that we are interested in finding connected branchs.
///  List Of  Residues:  All the residue nums of the branching from children residues
///  Tips:  All 'ends' of all the branches found using this function.
///
///  See Also: get_carbohydrate_residues_and_tips_of_branch
///            trim_carbohydrate_branch_from_X
void
get_branching_residues(
	Pose const & pose,
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
///  See Also: get_carbohydrate_residues_and_tips_of_branch
///            trim_carbohydrate_branch_from_X
void
fill_downstream_children_res_and_tips(
	Pose const & pose,
	Size res,
	Size parent_residue,
	utility::vector1< Size > & children_residues,
	utility::vector1< Size > & list_of_residues,
	utility::vector1< Size > & tips );

/// @brief Get all residue numbers in order from the tip to (and not including) stop_at_residue or a branch point.
///  All residue numbers are the tip or a linear polymer of glycans.
///  Useful for glycan stripping.
utility::vector1< core::Size >
get_resnums_in_leaf( Pose const & pose, Size tip_residue, Size stop_at_residue );

/// @brief Get all residue numbers in order from the tip to (and not including) stop_at_residue or a branch point.
///  All residue numbers are the tip or a linear polymer of glycans.
///  Useful for glycan stripping.
/// Do not use glycan_tree_set to get parent information.
utility::vector1< core::Size >
get_resnums_in_leaf_on_the_fly( Pose const & pose, Size tip_residue, Size stop_at_residue );

/// @brief Delete a leaf of glycan residues.  Use the ReferencePose to associate residue numbers.
/// @details
///  Uses delete_residue_slow as the debug_assert in delete_polymer_residue causes a crash.
///  This is a bug in the FoldTree where chemical edges are being treated as Jumps.
///  This is being addressed in a branch and a FoldTree.
void
delete_leaf( Pose & pose, utility::vector1< Size > leaf, std::string ref_pose_name = "temp_ref_pose" );


/// @brief Remove All/Any branch points from a carbohydrate or aa residue.
void
remove_carbohydrate_branch_point_variants( Pose & pose, core::Size const seqpos );

}  // namespace carbohydrates
}  // namespace pose
}  // namespace core

#endif  // INCLUDED_core_pose_carbohydrates_util_HH
