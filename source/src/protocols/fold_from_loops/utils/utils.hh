// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utils.hh
/// @brief  Helper functions for nub initio
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

#ifndef INCLUDED_protocols_fold_from_loops_utils_utils_hh
#define INCLUDED_protocols_fold_from_loops_utils_utils_hh

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

#include <string>

namespace protocols {
namespace fold_from_loops {
namespace utils {

/// @brief It will take each range that is not the first or the last one
/// and split it into two, following the logic of ```find_cutpoint_from_secondary_structure```
/// Input Ranges are not modified, a new ResidueRanges is returned.
/// WARNING1: When splitting a range of 1 residue, it splits into 1 range of 1 and another from 0 to 0, so
/// that it is easy to filter afterwards ( for the purpose of NubInitio, this means adding and empty pose
/// to the vector, for example ).
/// WARNING2: If the first ResidueRange does not start with 1, an empty_range (0,0) is set at the begining too.
core::select::residue_selector::ResidueRanges split_mid_ranges(
	std::string structure,
	core::select::residue_selector::ResidueRanges const & ranges );

/// @brief Given a Secondary Structure definition as a string, it finds what
/// it considers the optimal cutpoint.
/// This means that, if there is no loop 'L' segment it will return the middle
/// point as a residue; otherwise it will try to split the larger loop 'L' region.
core::Size find_cutpoint_from_secondary_structure( std::string structure );

/// @brief To a given scaffold it attaches a unfolded version of a pose to the N-terminal and another to the C-terminal
/// Generates a FoldTree that goes from the middle of 'scaffold' towards both directions
void attach_n_and_c_unfolded_poses_to_pose( core::pose::Pose const & n_insert, core::pose::Pose & scaffold, core::pose::Pose const & c_insert );
/// @brief adds the 'insert' unfolded pose to the N-terminal of the provided pose
void attach_unfolded_pose_to_pose_n_term( core::pose::Pose const & insert, core::pose::Pose & scaffold, core::chemical::ResidueTypeSetCOP rsd_set );
/// @brief adds the 'insert' unfolded pose to the C-terminal of the provided pose
void attach_unfolded_pose_to_pose_c_term( core::pose::Pose const & insert, core::pose::Pose & scaffold, core::chemical::ResidueTypeSetCOP rsd_set );
/// @brief Basically works as 'core::pose::append_pose_to_pose' but the FoldTree is build by joining by jump the
/// root of 'insert' FoldTree to that of 'scaffold' FoldTree. PDBInfo is lost.
void append_pose_to_pose_keep_fold_tree( core::pose::Pose & scaffold, core::pose::Pose const & insert, bool new_chain );

/// @brief Creates a sequence mapping betwee two proteins assuming that the False selections on a ResidueSubset
/// marks the residues that are the same for both of them. Assumes the same number of non-selected patches for each pose.
core::id::SequenceMapping map_by_residue_subsets( core::pose::Pose const & p1, core::select::residue_selector::ResidueSubset const & r1,
	core::pose::Pose const & p2, core::select::residue_selector::ResidueSubset const & r2);

/// @ brief reports the conditions of the pose
void report_unfolded( core::pose::Pose const & pose, core::kinematics::MoveMapOP movemap );
}
}
}

#endif
