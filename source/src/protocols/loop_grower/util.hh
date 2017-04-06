// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief utilities for Loop Growing Protocol
/// @details
/// @author Brandon Frenz


#ifndef INCLUDED_protocols_loop_grower_util_hh
#define INCLUDED_protocols_loop_grower_util_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.hh>

//// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace loop_grower {

//transforms an individual residue to all symmetric copies
void
transform_res_to_subunit( core::pose::Pose &pose, core::conformation::ResidueOP xformed, core::Size symmcopy);

//finds the closest symmetric unit and transforms all residues in sequence space to that position
void
transform_to_closest_symmunit(core::pose::Pose & cen_pose, core::pose::Pose & fa_pose, core::Size lower_pose);
}
}

#endif
