// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/grafting/util.hh
/// @brief Header for grafting utility functions.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Steven Lewis (smlewi@gmail.com)


#ifndef INCLUDED_protocols_grafting_util_hh
#define INCLUDED_protocols_grafting_util_hh

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

#include <protocols/loops/Loops.hh>

namespace protocols {
namespace grafting {


/// @brief Deletes a region of the pose. Starting from and including 'start' and 'end' residue.
void
delete_region(core::pose::Pose & pose, core::Size const start, core::Size const end);

/// @brief Returns a region of the pose including start and end as a new pose. Has a simple foldtree.
core::pose::Pose
return_region(core::pose::Pose & pose, core::Size const start, core::Size const end);

/// @brief replaces residues from from_pose to to_pose into pose where insertion region is defined. Returns product as a new value.
core::pose::Pose
replace_region(
	core::pose::Pose const & from_pose,
	core::pose::Pose const & to_pose,
	core::Size const from_pose_start_residue,
	core::Size const to_pose_start_residue,
	core::Size const insertion_length,
	bool copy_pdbinfo = false);

/// @author Steven Lewis smlewi@gmail.com, Jared Adolf-Bryfogle
/// @brief inserts one pose into another pose, returning the product as a new value.
/// @details Nter->Cter. Coordinates and dihedrals of insert are unchanged.
/// @details Begins insertion AFTER insert point.
core::pose::Pose
insert_pose_into_pose(
	core::pose::Pose const & scaffold_pose,
	core::pose::Pose const & insert_pose,
	core::Size const insert_point,
	core::Size const insert_point_end,
	bool copy_pdbinfo = false);

/// @brief inserts one pose into another pose, returning the product as a new value.
/// @details Nter->Cter. Coordinates and dihedrals of insert are unchanged.
/// @details Begins insertion AFTER insert point. insert_point_end is assumed to be insert_point+1.
core::pose::Pose
insert_pose_into_pose(
	core::pose::Pose const & scaffold_pose,
	core::pose::Pose const & insert_pose,
	core::Size const insert_point,
	bool copy_pdbinfo = false);


////////////////////////////////


/// @brief convenience function for AFTER apply method.
/// @details flexible Nter and Cter residues plus the first and last residue of the insert.
void
repack_connection_and_residues_in_movemap(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP fa_scorefxn,
	core::Size const start,
	core::Size const end,
	core::kinematics::MoveMapCOP movemap);

/// @brief convenience function for AFTER apply method.
/// @details flexible Nter and Cter residues plus the entire insert.
void
repack_connection_and_residues_in_movemap_and_piece(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP fa_scorefxn,
	core::Size const start,
	core::Size const end,
	core::kinematics::MoveMapCOP movemap);

/// @brief convenience function for AFTER apply method.
/// @details flexible Nter and Cter residues plus the entire insert and neighbors.
void
repack_connection_and_residues_in_movemap_and_piece_and_neighbors(
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionCOP fa_scorefxn,
	core::Size const start,
	core::Size const end,
	core::kinematics::MoveMapCOP movemap,
	core::Real neighbor_dis = 4.0);

/// @brief uses rms_util to superimpose overhang residues of piece onto pose.
/// @details Start + End denote residue number before and after the insert will be.
/// For example, start = 10, end = 11 for a scaffold where the previous residues are already deleted
/// or a scaffold where you are superposimposing a linker between two domains -
/// one that ends at start and the other that begins at end
void
superimpose_overhangs_heavy(
	core::pose::Pose const & pose,
	core::pose::Pose & piece,
	bool ca_only,
	core::Size start,
	core::Size end,
	core::Size Nter_overhang_len=2,
	core::Size Cter_overhang_length_len=2);

/// @brief deletes overhang residues of the pose piece.
/// Recommended use is within apply method
void
delete_overhang_residues(
	core::pose::Pose & piece,
	core::Size Nter_overhang_len,
	core::Size Cter_overhang_length_len);


/// @brief combines the two main movemaps to use after the insertion.
/// @details Start + End denote residue number before and after the insert.
/// original_end denotes the end residue number before insertion occurred
core::kinematics::MoveMapOP
combine_movemaps_post_insertion(
	core::kinematics::MoveMapCOP scaffold_mm,
	core::kinematics::MoveMapCOP insert_mm,
	core::Size start,
	core::Size original_end,
	core::Size insertion_length,
	core::Size cter_overhang_start = 0);

/// @brief Uses a small mover at high KT to perturb residues in the movemap for testing.
///  Returns bb_RMS_including_o
core::Real
perturb_backbone_for_test(core::pose::Pose & pose, core::kinematics::MoveMapOP mm);

/// @brief Idealize loop residues and residues in movemap.
// Optionally idealize insert residues as well.
void
idealize_combined_pose(
	core::pose::Pose & combined,
	core::kinematics::MoveMapOP movemap,
	core::Size start,
	core::Size insert_start,
	core::Size insert_end,
	core::Size Nter_loop_start,
	core::Size Cter_loop_end,
	bool idealize_insert = false);

/// @brief Adds cutpoint varients above and below cutpoint
void
add_cutpoint_variants_for_ccd(
	core::pose::Pose & pose,
	protocols::loops::Loops const & loops);

/// @brief Removes cutpoint variants above and below cutpoint
void
remove_cutpoint_variants_for_ccd(
	core::pose::Pose & pose,
	protocols::loops::Loops const & loops);

/// @brief Uses has_severe_peptide_bond_issues with stringent geometry values to
/// determine graft closure at cutpoint.
bool
graft_closed(
	core::pose::Pose & pose,
	protocols::loops::Loops & loops);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Useful FoldTrees
//
//

/////////////////////////////////////////////////////////////////
/// @brief ****Nter_loop_start---->Piece----> | Cter_loop_end****
///  Default FoldTree used by AnchoredGraftMover.
/// @params lower_cutpoint for CCD and loops is Cter_loop_end-1
///
void
setup_single_loop_single_arm_remodeling_foldtree(
	core::pose::Pose & pose,
	core::Size const Nter_loop_start,
	core::Size const Cter_loop_end,
	bool loop_modeling=false);

//////////////////////////////////////////////////////////////////
/// @brief ****Nter_loop_start---->Piece | <----Nter_loop_end****
/// Insert will move in cartesian space
/// @params lower_cutpoint for CCD and loops is end_-1
///
void
setup_single_loop_double_arm_remodeling_foldtree(
	core::pose::Pose & pose,
	core::Size const Nter_loop_start,
	core::Size const Cter_loop_end,
	core::Size end, bool loop_modeling=false);


}//namespace grafting
}//namespace protocols

#endif //INCLUDED_protocols_grafting_util_hh

