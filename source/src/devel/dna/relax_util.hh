// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Functions to help setup DNA relax.
/// 
/// @brief

#ifndef INCLUDED_devel_dna_relax_util_hh
#define INCLUDED_devel_dna_relax_util_hh


#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/kinematics/FoldTree.fwd.hh>
#include <utility/vector1.hh>


namespace devel {
namespace dna {


/// @brief  Add constraints to the pose's constraint set that try to keep the dna chain connected
void
setup_dna_chainbreak_constraints( core::pose::Pose & pose );


/// @brief  Choose a random base pair, returns the sequence number of the base-partner that comes first
core::Size
choose_random_base_pair( core::pose::Pose const & pose );


/// @brief  Choose a random base step jump, returns seqpos of 1st residue, not the jump number
core::Size
choose_random_base_step_jump( core::pose::Pose const & pose );


/// @brief  Adds jumps into the foldtree to support base-centric kinematics
void
add_dna_base_jumps_to_fold_tree(
																core::pose::Pose const & pose,
																core::kinematics::FoldTree & f,
																bool const flip=false
																);


/// @brief  Sets the jump-atoms in the foldtree so that the atomtree will have desired connectivity for dna jumps.
void
set_dna_jump_atoms( core::pose::Pose & pose );


/// @brief  Sets a foldtree for base-centric kinematics in a pose (legacy code)
void
setup_dna_only_fold_tree( core::pose::Pose & jump_pose, bool const flip = false );


/// @brief  Sets up a DNA-only pose by taking the paired residues in the first DNA chain in start_pose,
/// @brief  together with their partners. Sets the intra-base jumps in the foldtree to support base-centric kinematics.
void
setup_dna_only_jump_pose( core::pose::Pose const & start_pose, core::pose::Pose & jump_pose );


/// @brief  Delete unpaired DNA bases from the input pose
void
delete_unpaired_bases( core::pose::Pose & pose );

} // ns dna
} // ns devel

#endif
