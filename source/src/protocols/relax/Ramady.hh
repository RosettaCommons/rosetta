// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/relax/Ramady.hh
/// @brief Header for the Rana energy repair code, Ramady
/// @author Mike Tyka

#ifndef INCLUDED_protocols_relax_Ramady_hh
#define INCLUDED_protocols_relax_Ramady_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace relax {

void add_coordinate_constraints_to_pose( core::pose::Pose & pose, const core::pose::Pose &constraint_target_pose,  protocols::loops::Loops &exclude_regions );
void fix_worst_bad_ramas( core::pose::Pose & original_pose, core::Size how_many = 1, core::Real skip_prob = 0.0, core::Real limit_rms=0.5, core::Real limit_rama = 2.0 );


} // relax
} // protocols

#endif
