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
/// @author


#ifndef INCLUDED_protocols_loops_loop_closure_ccd_ccd_closure_hh
#define INCLUDED_protocols_loops_loop_closure_ccd_ccd_closure_hh


// Package headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

int
fast_ccd_loop_closure(
	core::pose::Pose & pose,
	core::kinematics::MoveMap const & mm,
	int const loop_begin,
	int const loop_end,
	int const cutpoint,
	int const ii_cycles,
	core::Real const tolerance,
	bool const rama_check_boltzmann,
	core::Real const max_rama_score_increase,
	core::Real const max_total_delta_helix,
	core::Real const max_total_delta_strand,
	core::Real const max_total_delta_loop,
	core::Real & forward_deviation, // output
	core::Real & backward_deviation, // output
	core::Real & torsion_delta,
	core::Real & rama_delta
);

void
ccd_moves(
	int const total_moves,
	core::pose::Pose & pose,
	core::kinematics::MoveMap const & mm,
	int const loop_begin,
	int const loop_end,
	int const cutpoint
);

std::pair<core::Real, core::Real>
get_deviation(
	core::pose::Pose const & pose,
	int const cutpoint
);

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_loop_closure_ccd_ccd_closure_hh
