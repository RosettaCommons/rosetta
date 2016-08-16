// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @details
/// @author


#ifndef INCLUDED_protocols_loops_loopfinder_HH
#define INCLUDED_protocols_loops_loopfinder_HH

// Package headers
#include <protocols/loops/Loops.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace loops {

//@brief finds protein loops from dssp definitions - useful if you want to redesign all loops in a pose
void loopfinder(  core::pose::Pose & pose , loops::Loops & loops );

} //loops
} //protocols

#endif //INCLUDED_protocols_loops_loopfinder_HH
