// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/movemap/util.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_movemap_util_HH
#define INCLUDED_protocols_stepwise_sampling_movemap_util_HH

#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/toolbox/AllowInsert.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace sampling {

void
figure_out_stepwise_movemap( core::kinematics::MoveMap & mm,
														 core::pose::Pose const & pose,
														 utility::vector1< core::Size > const & working_minimize_res,
														 bool const move_takeoff_torsions = true );

void
figure_out_stepwise_movemap( core::kinematics::MoveMap & mm,
														 toolbox::AllowInsertOP & allow_insert,
														 core::pose::Pose const & pose,
														 utility::vector1< core::Size > const & working_fixed_res,
														 utility::vector1< core::Size > const & working_extra_minimize_res,
														 bool const move_takeoff_torsions = true );

void
figure_out_stepwise_movemap( core::kinematics::MoveMap & mm,
														 toolbox::AllowInsertOP allow_insert,
														 core::pose::Pose const & pose,
														 utility::vector1< core::Size > const & working_minimize_res,
														 bool const move_takeoff_torsions = true );

	// deprecated?
void
figure_out_moving_residues( core::kinematics::MoveMap & mm, core::pose::Pose const & pose,
														utility::vector1< core::Size > const & fixed_res,
														bool const move_takeoff_torsions = true,
														bool const move_jumps_between_chains = false );

void
check_move_map_against_working_parameters( core::pose::Pose const & pose ,
																					 core::kinematics::MoveMapCOP minimize_map,
																					 working_parameters::StepWiseWorkingParametersCOP working_parameters );



} //sampling
} //stepwise
} //protocols

#endif
