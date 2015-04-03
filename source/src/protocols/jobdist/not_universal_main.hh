// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jobdist/not_universal_main.hh
/// @brief  Simple main method for applying a Mover to a set of
/// input Poses.
/// @author James Thompson

#ifndef INCLUDED_protocols_jobdist_not_universal_main_hh
#define INCLUDED_protocols_jobdist_not_universal_main_hh

#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// Utility headers

// C++ headers

#include <utility/vector1.hh>


namespace protocols {
namespace jobdist {

bool pose_matches_user_tag(
	core::pose::Pose & pose,
	utility::vector1< std::string > const & user_tags
);

int not_universal_main(
	protocols::moves::Mover & mover
);

} // namespace jobdist
} // namespace protocols

#endif // INCLUDED_protocols_jobdist_not_universal_main_HH
