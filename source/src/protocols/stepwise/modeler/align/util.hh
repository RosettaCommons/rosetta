// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/align/util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_align_util_HH
#define INCLUDED_protocols_stepwise_modeler_align_util_HH

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <utility/vector1.fwd.hh>
#include <map>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace align {

core::Real
get_rmsd( core::pose::Pose const & pose1, core::pose::Pose const & pose2,
	utility::vector1< core::Size > const & calc_rms_res,
	bool const check_align_at_superimpose_res = false,
	bool const check_switch = false );

core::Real
get_rmsd( core::pose::Pose const & pose1, core::pose::Pose const & pose2,
	bool const check_align_at_superimpose_res = false,
	bool const check_switch = false );


void
align_pose_and_add_rmsd_constraints( core::pose::Pose & pose,
	core::pose::PoseCOP align_pose,
	utility::vector1< core::Size > const & moving_res_list,
	core::Real const rmsd_screen );

core::Real
superimpose_with_stepwise_aligner( core::pose::Pose & pose, core::pose::Pose const & align_pose,
	bool superimpose_over_all_instantiated = false );

///////////////////////////////////////////////////////////////////////////////////////////
// Following functions (superimpose_pose, creat_alignment_id_map) use legacy code
// for choosing which atoms to superimpose on -- but are called by InputStreamWithResidueInfo
// and a couple other classes that should probably ALL BE DEPRECATED. -- rhiju, 2014
///////////////////////////////////////////////////////////////////////////////////////////
/// @brief  Superimpose mod_pose onto ref_pose using the mapping of residues from
/// mod_pose to ref_pose given by res_map
core::Real
superimpose_pose_legacy(
	core::pose::Pose & mod_pose,
	core::pose::Pose const & ref_pose,
	std::map< core::Size, core::Size > const & res_map
);

// This should be deprecated in favor of StepWisePoseAligner
core::id::AtomID_Map< core::id::AtomID >
create_alignment_id_map_legacy( core::pose::Pose const & mod_pose, core::pose::Pose const & ref_pose,
	utility::vector1< core::Size > const & superimpose_res );

// This should be deprecated in favor of StepWisePoseAligner
core::id::AtomID_Map< core::id::AtomID >
create_alignment_id_map_legacy( core::pose::Pose const & mod_pose,
	core::pose::Pose const & ref_pose,
	std::map< core::Size, core::Size > res_map );


} //align
} //modeler
} //stepwise
} //protocols

#endif
