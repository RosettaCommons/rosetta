// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  Header file for the solution picker functions.
/// @author Kale Kundert (kale.kundert@ucsf.edu)

#ifndef INCLUDED_protocols_kinematic_closure_utilities_HH
#define INCLUDED_protocols_kinematic_closure_utilities_HH

// Unit headers
#include <protocols/kinematic_closure/ClosureProblem.fwd.hh>
#include <protocols/kinematic_closure/ClosureSolution.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>

namespace protocols {
namespace kinematic_closure {

/// @brief Setup the given pose with a fold tree that is ideally configured for 
/// sampling the given loop with kinematic closure.

void setup_fold_tree(
		core::pose::Pose & pose,
		protocols::loops::Loop const & loop);

}
}

#endif



