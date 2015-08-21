// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/sample_around/util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_toolbox_sample_around_util_HH
#define INCLUDED_protocols_toolbox_sample_around_util_HH

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>


namespace protocols {
namespace toolbox {
namespace sample_around {

void
add_virtual_res ( core::pose::Pose & pose, bool set_res_as_root = true );

void
add_another_virtual_res ( core::pose::Pose & pose );

void
rotate_into_nucleobase_frame( core::pose::Pose & pose );

void
rotate_into_phosphate_frame( core::pose::Pose & pose, core::Size const n, bool const center_on_OP2 );

core::Real
centroid_dist( core::pose::Pose & pose,
	bool const sample_another_adenosine = false );

core::Real
sample_all_rotations_at_jump( core::pose::Pose & pose, core::Size const num_jump, core::scoring::ScoreFunctionOP scorefxn = 0 );

core::Real
do_scoring( core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP scorefxn,
	bool const & sample_rotations,
	core::Size const probe_jump_num );

void
do_xy_scan( core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP scorefxn,
	std::string const & outfile,
	core::Real const z,
	core::Size const probe_jump_num,
	core::Real const box_bins,
	core::Real const translation_increment,
	bool const sample_rotations );

} //sample_around
} //toolbox
} //protocols

#endif
