// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/util.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_rna_StepWiseMonteCarloUtil_HH
#define INCLUDED_protocols_stepwise_monte_carlo_rna_StepWiseMonteCarloUtil_HH

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace monte_carlo {

void
output_to_silent_file( std::string const & out_tag,
											 std::string const & silent_file,
											 pose::Pose & pose,
											 pose::PoseCOP native_pose,
											 bool const superimpose_over_all_instantiated = false,
											 bool const do_rms_fill_calculation = false );

core::io::silent::SilentStructOP
prepare_silent_struct( std::string const & out_tag,
											 pose::Pose & pose,
											 pose::PoseCOP native_pose,
											 bool const superimpose_over_all_instantiated = false,
											 bool const do_rms_fill_calculation = false,
											 core::pose::PoseOP full_model_pose  = 0 );

void
output_to_silent_file( std::string const & out_tag,
											 std::string const & silent_file,
											 pose::Pose const & pose );

void
output_to_silent_file( std::string const & silent_file,
											 utility::vector1< pose::PoseOP > & pose_list,
											 pose::PoseCOP native_pose );

void
build_full_model( core::pose::Pose const & start_pose, core::pose::Pose & full_model_pose );

core::pose::PoseOP
build_full_model( pose::Pose const & start_pose );

} //monte_carlo
} //stepwise
} //protocols

#endif
