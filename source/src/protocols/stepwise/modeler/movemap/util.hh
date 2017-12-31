// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/movemap/util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_movemap_util_HH
#define INCLUDED_protocols_stepwise_modeler_movemap_util_HH

#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace movemap {

void
figure_out_stepwise_movemap( core::kinematics::MoveMap & mm,
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & working_minimize_res,
	bool const move_takeoff_torsions = true );

void
figure_out_stepwise_movemap( core::kinematics::MoveMap & mm,
	core::pose::toolbox::AtomLevelDomainMapOP & atom_level_domain_map,
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & working_fixed_res,
	utility::vector1< core::Size > const & working_extra_minimize_res,
	bool const move_takeoff_torsions = true );

void
figure_out_stepwise_movemap( core::kinematics::MoveMap & mm,
	core::pose::toolbox::AtomLevelDomainMapOP atom_level_domain_map,
	core::pose::Pose const & pose,
	utility::vector1< core::Size > const & working_minimize_res,
	bool const move_takeoff_torsions = true );

} //movemap
} //modeler
} //stepwise
} //protocols

#endif
