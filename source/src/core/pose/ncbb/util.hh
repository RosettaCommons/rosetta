// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/pose/ncbb/util.hh
/// @brief   Utility function declarations for poses with noncanonical backbone. 
/// @author  kdrew
/// @author  Andy Watkins

#ifndef INCLUDED_core_pose_ncbb_util_HH
#define INCLUDED_core_pose_ncbb_util_HH

// Unit header
#include <core/pose/Pose.fwd.hh>

// Project headers
#include <core/types.hh>
//#include <core/conformation/Residue.fwd.hh>
//#include <core/conformation/Conformation.fwd.hh>

// C++ headers
//#include <string>


namespace core {
namespace pose {
namespace ncbb {

// @brief initializes ncbbs in pose, returns residue numbers of oops
utility::vector1< core::Size > initialize_ncbbs( Pose & pose); 
// @brief initializes oops in pose, returns residue numbers of oops
utility::vector1< core::Size > initialize_oops( Pose & pose); 
utility::vector1< core::Size > initialize_hbs( Pose & pose); 
/// @brief  Add constraints to keep oligooxopiperazine (oop) ring closed, default values (distance = 1.5, std = 0.05)
void add_oop_constraint( core::pose::Pose & pose, core::Size oop_seq_position );
/// @brief  Add constraints to keep oligooxopiperazine (oop) ring closed 
void add_oop_constraint( core::pose::Pose & pose, core::Size oop_seq_position, core::Real distance, core::Real std );
/// @brief  Add constraints to keep hydrogen bond surrogate (hbs) macrocycle closed, default values (distance = 1.52, std = 0.05)
void add_hbs_constraint( core::pose::Pose & pose, core::Size oop_seq_position );
/// @brief  Add constraints to keep hydrogen bond surrogate (hbs) ring closed 
void add_hbs_constraint( core::pose::Pose & pose, core::Size hbs_seq_position, core::Real distance, core::Real std );

}  // namespace ncbb
}  // namespace pose
}  // namespace core

#endif  // INCLUDED_core_pose_ncbb_util_HH
