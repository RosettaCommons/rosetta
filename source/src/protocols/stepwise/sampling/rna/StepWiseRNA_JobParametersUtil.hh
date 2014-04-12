// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/StepWiseRNA_JobParametersUtil.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_rna_StepWiseRNA_JobParametersUtil_HH
#define INCLUDED_protocols_stepwise_sampling_rna_StepWiseRNA_JobParametersUtil_HH

#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParametersUtil.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/toolbox/AllowInsert.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <utility/vector1.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {

StepWiseRNA_JobParametersOP
setup_job_parameters_for_swa( Size const moving_res,
															pose::Pose const & pose,
															pose::PoseCOP native_pose,
															utility::vector1< Size > const & calc_rms_res_,
															utility::vector1< Size > const & terminal_res_,
															utility::vector1< Size > const & syn_chi_res_list_,
															utility::vector1< Size > const & extra_minimize_res_,
															utility::vector1< Size > & fixed_res_,
															utility::vector1< Size > & minimize_res_,
															toolbox::AllowInsertOP & allow_insert_,
															kinematics::MoveMapOP & minimize_move_map );

void
figure_out_root_partition_res( pose::Pose const & pose,
															 StepWiseRNA_JobParametersCOP job_parameters,
															 utility::vector1< Size > & root_partition_res,
															 utility::vector1< Size > & moving_partition_res );

StepWiseRNA_JobParametersOP
setup_job_parameters_legacy( Size const rebuild_res,
														 pose::Pose const & pose,
														 pose::PoseCOP native_pose,
														 utility::vector1< Size > const & calc_rms_res,
														 utility::vector1< Size > const & terminal_res,
														 utility::vector1< Size > const & syn_chi_res_list,
														 utility::vector1< Size > & fixed_res_guess,
														 utility::vector1< Size > & suites_that_must_be_minimized );

StepWiseRNA_JobParametersOP
setup_job_parameters_explicit( Size const rebuild_res,
															 pose::Pose const & pose,
															 pose::PoseCOP native_pose,
															 utility::vector1< Size > const & calc_rms_res,
															 utility::vector1< Size > const & terminal_res,
															 utility::vector1< Size > const & syn_chi_res_list,
															 utility::vector1< Size > & fixed_res_guess,
															 utility::vector1< Size > & suites_that_must_be_minimized );

bool
figure_out_rebuild_bulge_mode( pose::Pose const & pose, Size const rebuild_res );

bool
figure_out_sample_both_sugar_base_rotamer( pose::Pose const & pose, bool const floating_base, Size const rebuild_suite );

// Undefined, commenting out to fix PyRosetta build  bool
// Undefined, commenting out to fix PyRosetta build  figure_out_is_residue_prepend( Size const seq_num, StepWiseRNA_JobParametersOP job_parameters );

} //rna
} //sampling
} //stepwise
} //protocols

#endif
