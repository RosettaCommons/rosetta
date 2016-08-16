// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/kinematics/visualize.hh
/// @brief  3-D visualizations of FoldTree and AtomTree in kinemage format
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_viewer_visualize_hh
#define INCLUDED_protocols_viewer_visualize_hh


// Package headers

// Rosetta headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// ObjexxFCL headers

// C++ Headers
#ifdef WIN32
#include <string>
#endif

namespace protocols {
namespace viewer {


void
dump_pose_kinemage(
	std::string const & filename,
	core::pose::Pose const & pose
);


} // namespace viewer
} // namespace protocols


#endif // INCLUDED_core_kinematics_visualize_HH
