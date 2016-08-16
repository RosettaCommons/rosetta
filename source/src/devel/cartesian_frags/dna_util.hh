// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#ifndef INCLUDED_devel_cartesian_frags_dna_util_hh
#define INCLUDED_devel_cartesian_frags_dna_util_hh


// libRosetta headers

#include <devel/cartesian_frags/CartesianFragment.fwd.hh>
#include <devel/cartesian_frags/DNA_FragLib.fwd.hh>


#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>
#include <devel/cartesian_frags/Direction.hh>
#include <utility/vector1.hh>


namespace devel {
namespace cartesian_frags {

/// @brief  Optimize the bond angles and torsions of a sugar->sugar backbone suite to match a desired transform

void
optimize_suite(
	Direction const & dir,
	core::kinematics::RT const & target,
	core::pose::Pose & pose, // the mini pose
	CartesianFragment & frag
);


/// @brief  Scan a fragment library for fragment combinations that patch up the DNA backbone.

core::Real
find_sugar_and_suite_frags(
	DNA_FragLib const & lib,
	core::kinematics::RT const & fwd_target,
	core::kinematics::RT const & bwd_target,
	bool const lower_terminus,
	bool const upper_terminus,
	CartesianFragment & sugar,
	CartesianFragment & fwd_suite,
	CartesianFragment & bwd_suite
);


/// @brief  Patches up the backbone before and/or after DNA position i by varying the sugar and both suites

core::Real
patch_up_backbone(
	core::Size const i,
	DNA_FragLib const & lib,
	core::conformation::Conformation & conf
);

/// @brief  Patches up the DNA backbone link between i and i+1
void
patch_up_backbone_link(
	core::Size const i,
	DNA_FragLib const & lib,
	core::scoring::ScoreFunction const & scorefxn,
	core::pose::Pose & pose_inout
);


} // ns cartesian_frags
} // ns devel

#endif
