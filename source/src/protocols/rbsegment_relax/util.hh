// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Frank DiMaio
/// @author Srivatsan Raman

#ifndef INCLUDED_protocols_rbsegment_relax_util_hh
#define INCLUDED_protocols_rbsegment_relax_util_hh

// Package headers
#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <protocols/rbsegment_relax/RBSegment.fwd.hh>
#include <protocols/rbsegment_relax/RBSegment.hh>

#include <map>

#include <numeric/xyzVector.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <numeric/conversions.hh>
#include <numeric/xyzMatrix.hh>


namespace protocols {
namespace rbsegment_relax {

/// @brief set up constraints over RB segments only; allow ambiguity in sequence threading
void set_rb_constraints(
	core::pose::Pose & pose,
	core::pose::Pose const &cst_pose,
	utility::vector1< protocols::rbsegment_relax::RBSegment > const & rbsegs ,
	core::id::SequenceMapping const & resmap,
	core::Real cst_width,
	core::Real cst_stdev,
	core::Size cst_seqwidth );

/// @brief set up constraints accounting for missing density in start pose
void set_constraints(
	core::pose::Pose & pose,
	core::pose::Pose const &cst_pose,
	core::Real cst_width,
	core::Real cst_stdev,
	core::Size cst_seqwidth );

/// @brief setup star-topology fold tree
void setup_star_topology( core::pose::Pose & pose );

/// @brief build a pose where every residue is connected by jumps to vrt
void setup_disconnected( core::pose::Pose & pose );

/// @brief remove loops from pose
void setup_pose_from_rbsegs(
	utility::vector1< protocols::rbsegment_relax::RBSegment > const &rbsegs ,
	core::pose::Pose const &pose_in ,
	core::pose::Pose &pose_out ,
	core::id::SequenceMapping &resmap,
	core::kinematics::MoveMap &mm ,
	bool fixligs=false );


/// @brief use DSSP and simple rules to guess the asignment of rigid-body segments
void guess_rbsegs_from_pose(
	core::pose::Pose const & pose,
	utility::vector1< RBSegment > & rigid_segs,
	utility::vector1< RBSegment > & rb_chunks,
	protocols::loops::Loops & loops
);

/// @brief
utility::vector1<core::Size> setup_pose_rbsegs_keep_loops(
	core::pose::Pose &pose,
	utility::vector1< protocols::rbsegment_relax::RBSegment > const &rbsegs ,
	protocols::loops::Loops const &loops,
	core::kinematics::MoveMapOP mm );

/// @brief restore loops from pose
void restore_pose_from_rbsegs(
	utility::vector1< protocols::rbsegment_relax::RBSegment > const &rbsegs ,
	core::pose::Pose const &pose_in ,
	core::pose::Pose &pose_out /* input/output */ );

/// @apply res mapping to rbsegments
void remap_rb_segments(
	utility::vector1< RBSegment > const &rbsegs,
	utility::vector1< RBSegment > &rbsegs_remap,
	core::id::SequenceMapping const &resmap);


}
}

#endif
