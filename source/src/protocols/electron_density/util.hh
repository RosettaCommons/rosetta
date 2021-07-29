// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief protocols for folding into density
/// @details
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_electron_density_util_hh
#define INCLUDED_protocols_electron_density_util_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/electron_density/ElectronDensity.fwd.hh>

//// C++ headers
#include <string>

#include <protocols/loops/Loops.fwd.hh>


namespace protocols {
namespace electron_density {

// docks the pose into the map using the protocol specified in -edensity::realign
core::Real dockPoseIntoMap( core::pose::Pose & pose , std::string const & align_in="");

// find N residues with worst agreement to density
protocols::loops::Loops findLoopFromDensity( core::pose::Pose & pose, core::Real frac, int max_helix, int max_strand );

void
remove_occupied_density_from_density(
	core::pose::Pose const & pose,
	core::scoring::electron_density::ElectronDensity & dens,
	core::Size const edge_trim = 5,
	core::Real const mask_radius = 2.0);

}
}

#endif
