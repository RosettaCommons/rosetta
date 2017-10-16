// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/variant_util.hh
/// @brief  Pose utilities
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov, Vikram K. Mulligan, Jared Adolf-Bryfogle

#ifndef INCLUDED_core_pose_variant_util_hh
#define INCLUDED_core_pose_variant_util_hh

// Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/MiniPose.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/rings/AxEqDesignation.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/io/StructFileRep.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/scoring/ScoreType.hh>

// Utility headers
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <map>
#include <set>

namespace core {
namespace pose {

/// @brief Remove variant from an existing residue.
conformation::ResidueOP remove_variant_type_from_residue(
	conformation::Residue const & old_rsd,
	core::chemical::VariantType const variant_type,
	pose::Pose const & pose );

/// @brief Construct a variant of an existing residue.
conformation::ResidueOP add_variant_type_to_residue(
	conformation::Residue const & old_rsd,
	core::chemical::VariantType const variant_type,
	pose::Pose const & pose );

/// @brief Construct a variant of an existing pose residue.
void add_variant_type_to_pose_residue(
	pose::Pose & pose,
	chemical::VariantType const variant_type,
	Size const seqpos );

/// @brief Construct a non-variant of an existing pose residue.
void remove_variant_type_from_pose_residue(
	pose::Pose & pose,
	chemical::VariantType const variant_type,
	Size const seqpos );

void
add_lower_terminus_type_to_pose_residue(
	pose::Pose & pose,
	Size const seqpos
);

void
add_upper_terminus_type_to_pose_residue(
	pose::Pose & pose,
	Size const seqpos
);

void
remove_lower_terminus_type_from_pose_residue(
	pose::Pose & pose,
	Size const seqpos
);

void
remove_upper_terminus_type_from_pose_residue(
	pose::Pose & pose,
	Size const seqpos
);

/// @brief Add cutpoint variants to all residues annotated as cutpoints in the FoldTree in the Pose.
void
correctly_add_cutpoint_variants( core::pose::Pose & pose );

/// @brief Add CUTPOINT_LOWER and CUTPOINT_UPPER types to two residues, remove incompatible types, and declare
/// a chemical bond between them.
/// @param[in] pose The pose to modify.
/// @param[in] cutpoint_res The index of the CUTPOINT_LOWER residue.
/// @param[in] check_fold_tree If true, a check is performed to confirm that the residues in question represent a
/// cutpoint in the foldtree in the pose.
/// @param[in] next_res_in The index of the CUTPOINT_UPPER residue.  If not provided, or if set to 0, this defaults
/// to the cutpoint_res + 1 residue.  Must be specified for cyclic geometry.
void
correctly_add_cutpoint_variants(
	core::pose::Pose & pose,
	Size const cutpoint_res,
	bool const check_fold_tree = true,
	Size const next_res_in = 0 );

/// @brief Remove variant types incompatible with CUTPOINT_LOWER from a position in a pose.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
/// @param[in,out] pose The pose on which to operate.
/// @param[in] res_index The index of the residue on which to operate.
void
correctly_remove_variants_incompatible_with_lower_cutpoint_variant( core::pose::Pose & pose, Size const res_index );

/// @brief Remove variant types incompatible with CUTPOINT_UPPER from a position in a pose.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
/// @param[in,out] pose The pose on which to operate.
/// @param[in] res_index The index of the residue on which to operate.
void
correctly_remove_variants_incompatible_with_upper_cutpoint_variant( core::pose::Pose & pose, Size const res_index );

/// @brief returns true if the given residue in the pose is a chain ending or has upper/lower terminal variants
bool
pose_residue_is_terminal( Pose const & pose, Size const resid );

/// @brief checks to see if this is a lower chain ending more intelligently than just checking residue variants
bool
is_lower_terminus( pose::Pose const & pose, Size const resid );

/// @brief checks to see if this is a lower chain ending more intelligently than just checking residue variants
bool
is_upper_terminus( pose::Pose const & pose, Size const resid );

void
fix_up_residue_type_variants_at_strand_beginning( core::pose::Pose & pose, core::Size const res );

void
fix_up_residue_type_variants_at_strand_end( core::pose::Pose & pose, core::Size const res );

void
fix_up_residue_type_variants( core::pose::Pose & pose );


} // pose
} // core

#endif // INCLUDED_core_pose_variant_util_HH
