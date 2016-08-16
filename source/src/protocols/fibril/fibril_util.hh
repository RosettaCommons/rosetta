// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/fibril/fibril_util.hh
/// @brief utility functions for handling with symmetric fibril
/// @author Lin Jiang

#ifndef INCLUDED_protocols_fibril_fibril_util_hh
#define INCLUDED_protocols_fibril_fibril_util_hh


// Unit headers


#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace fibril {

void
reorient_extended_fibril(
	core::conformation::Conformation & src_conformation,
	core::conformation::symmetry::SymmData & symmdata
);

void
make_symmetric_fibril(
	core::pose::Pose & pose
);

void
superimpose_pose_on_subset_bb(
	core::pose::Pose& pose,
	core::pose::Pose& ref_pose,
	protocols::loops::Loops core,
	protocols::loops::Loops ref_core
);


} // fibril
} // protocols


#endif
