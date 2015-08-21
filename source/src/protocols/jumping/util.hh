// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file JumpSetup
/// @brief read jump-definition file   setups fold tree an chainbreak variants
/// loop code didn't work because fold-tree to complicated ( overlapping loops )
/// @details
/// @author Oliver Lange


#ifndef INCLUDED_protocols_jumping_util_hh
#define INCLUDED_protocols_jumping_util_hh

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/loops/loop_closure/ccd/SlidingWindowLoopClosure.fwd.hh>
#include <core/types.hh>
#include <core/fragment/JumpingFrame.fwd.hh>

#include <protocols/checkpoint/CheckPointer.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace jumping {

//@brief remove cutpoint from fold-tree and idealize bond/angles at cutpoint and cutpoint+1
// note that pose_out will have weird structure --> use e.g., SlidingWindowLoopClosure to repair that
// returns false if there is a failure
bool remove_cut(
	core::Size cutpoint,
	core::pose::Pose& pose,
	core::kinematics::FoldTree const& final_fold_tree = core::kinematics::FoldTree()
);

bool remove_cut(
	core::Size cutpoint,
	core::kinematics::FoldTree &fold_tree,
	core::kinematics::FoldTree const& final_fold_tree = core::kinematics::FoldTree()
);

//@brief remove all cutpoints, close loops using the supplied closure_protocol
void close_chainbreaks(
	loops::loop_closure::ccd::SlidingWindowLoopClosureOP closure_protocol,
	core::pose::Pose& pose,
	checkpoint::CheckPointer &checkpoint,
	const std::string &decoy_tag,
	core::kinematics::FoldTree const& final_fold_tree = core::kinematics::FoldTree()
);

void safe_secstruct( core::pose::Pose& pose );


void
get_pleating(
	core::pose::Pose const& pose,
	core::Size const pos1,
	core::Size const pos2,
	core::Size &orientation,
	core::Size &pleating
);

void
get_pairing_geometry(
	core::pose::Pose const& pose,
	core::Size const res1,
	core::Size const res2,
	core::Real& orientation,
	core::Real& pleating1,
	core::Real& pleating2
);

//@brief genrates the standard SS-Jump Frame  BBTorsion UpJump DownJump BBTorsion
// depends on length if both BBTorsions are present...
core::fragment::JumpingFrameOP generate_empty_jump_frame(
	core::Size startpos,
	core::Size endpos,
	core::Size length
);

//@brief generate the standard SS-Jump Frame with a template FragData --- use this if you want to steal jumps
core::fragment::JumpingFrameOP generate_jump_frame(
	core::Size startpos,
	core::Size endpos,
	bool bWithTorsion
);

/// @brief Assign secondary structure using DSSP.
void assign_ss_dssp( core::pose::Pose & pose );

} //jumping
} //protocols

#endif
