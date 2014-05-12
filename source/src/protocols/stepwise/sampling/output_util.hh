// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/output_util.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_output_util_HH
#define INCLUDED_protocols_stepwise_sampling_output_util_HH

#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <iostream>

namespace protocols {
namespace stepwise {
namespace sampling {

void
output_boolean( std::string const & tag, bool boolean, std::ostream & outstream = std::cout );

void
output_boolean( bool boolean, std::ostream & outstream = std::cout );

void
output_movemap( core::kinematics::MoveMap const & mm, core::pose::Pose const & pose, std::ostream & outstream = std::cout );

} //sampling
} //stepwise
} //protocols

#endif
