// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/subpose_manipulation_util.hh
/// @brief  Pose subpose_manipulation_utilities
/// @author Phil Bradley
/// @author Modified by Sergey Lyskov, Vikram K. Mulligan, Jared Adolf-Bryfogle

#ifndef INCLUDED_core_pose_subpose_manipulation_util_hh
#define INCLUDED_core_pose_subpose_manipulation_util_hh

// Package headers
#include <core/pose/Pose.fwd.hh>
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

typedef std::set< int > Jumps;

/// @brief Append residues of pose2 to pose1.
void
append_pose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	bool new_chain = true
);

/// @brief Append specified residues of pose2 to pose1.
void
append_subpose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	core::Size start_res,
	core::Size end_res,
	bool new_chain = true
);

/// @brief Create a subpose of the src pose.  PDBInfo is set as NULL.
void
create_subpose(
	Pose const & src,
	utility::vector1< Size > const & positions,
	kinematics::FoldTree const & f,
	Pose & pose
);

/// @brief Create a subpose of the src pose -- figures out a reasonable fold tree.
void
pdbslice( pose::Pose & new_pose,
	pose::Pose const & pose,
	utility::vector1< Size > const & slice_res );

/// @brief Create a subpose of the src pose -- figures out a reasonable fold tree.
void
pdbslice( pose::Pose & pose,
	utility::vector1< Size > const & slice_res );

// for partition_by_jump: both new poses start residue numbering from 1 and don't keep the original numbering!
void
partition_pose_by_jump(
	pose::Pose const & src,
	int const jump_number,
	pose::Pose & partner1, // partner upstream in foldtree
	pose::Pose & partner2  // partner downstream in foldtree
);

} // pose
} // core

#endif // INCLUDED_core_pose_subpose_manipulation_util_HH
