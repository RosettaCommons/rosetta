// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/idealize/idealize.hh
/// @brief protocols for idealizing a Pose
/// @author


#ifndef INCLUDED_protocols_idealize_idealize_hh
#define INCLUDED_protocols_idealize_idealize_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1_bool.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace idealize {

void
dihedral_distance(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	utility::vector1< bool > const & use_pos,
	core::Real & avg_bb_angle_dev,
	core::Real & max_bb_angle_dev,
	core::Real & avg_chi_angle_dev,
	core::Real & max_chi_angle_dev
);

void
basic_idealize(
	core::pose::Pose & pose,
	utility::vector1< core::Size > pos_list, // local copy
	core::scoring::ScoreFunction const & scorefxn,
	bool const fast,
	bool const chainbreaks = false, // to maintain previous behavior (abrelax and ?  may use idealizer to close chainbreaks)
	bool const cis_omega = false //fpd fix non-proline cis-omegas
);

} // idealize
} // protocols

#endif
