// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief protocols for folding into density
/// @details
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_electron_density_util_hh
#define INCLUDED_protocols_electron_density_util_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

//// C++ headers
#include <string>

#include <protocols/loops/Loops.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace electron_density {

// docks the pose into the map using the protocol specified in -edensity::realign
core::Real dockPoseIntoMap( core::pose::Pose & pose , std::string align_in="");

// helper function called by dockPoseIntoMap
core::Real fastTransAlignPose( core::pose::Pose & pose );

// helper function called by dockPoseIntoMap -- 2D rotation alignment
core::Real fast2DRotAlignPose( core::pose::Pose & pose , std::string axis );

// find stretch of N residues with worst agreement to patterson
protocols::loops::Loops findLoopFromPatterson( core::pose::Pose & pose, core::Size N, core::Size nloops, bool allow_termini );

// find N residues with worst agreement to density
protocols::loops::Loops findLoopFromDensity( core::pose::Pose & pose, core::Real frac, int max_helix, int max_strand );

}
}

#endif
