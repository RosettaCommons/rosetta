// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/output_util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_output_util_HH
#define INCLUDED_protocols_stepwise_modeler_output_util_HH

#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <iostream>

namespace protocols {
namespace stepwise {
namespace modeler {

void
output_boolean( std::string const & tag, bool boolean, std::ostream & outstream = std::cout );

void
output_boolean( bool boolean, std::ostream & outstream = std::cout );

//following copies code that is in rna/StepWiseRNA_Util? Remove the latter?
void
output_boolean(std::string const & tag, bool boolean, std::ostream & TR );

void
output_boolean(bool boolean, std::ostream & TR );

void
output_movemap( core::kinematics::MoveMap const & mm, core::pose::Pose const & pose, std::ostream & outstream = std::cout );

} //modeler
} //stepwise
} //protocols

#endif
