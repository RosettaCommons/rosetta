// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_BetaAntiparallelUtil.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_protein_StepWiseBetaAntiparallelUtil_HH
#define INCLUDED_protocols_stepwise_protein_StepWiseBetaAntiparallelUtil_HH

#include <protocols/rotamer_sampler/RotamerSized.fwd.hh>
#include <protocols/stepwise/sampling/modeler_options/StepWiseModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/sampling/protein/InputStreamWithResidueInfo.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace rotamer_sampler {
namespace protein {

	RotamerSizedOP
	get_basic_protein_sampler(
      core::pose::Pose const & pose,
			utility::vector1< core::Size > const & moving_res_list,
			protocols::stepwise::sampling::working_parameters::StepWiseWorkingParametersCOP working_parameters,
			protocols::stepwise::sampling::modeler_options::StepWiseModelerOptionsCOP options,
			utility::vector1< protocols::stepwise::sampling::protein::InputStreamWithResidueInfoOP > & input_streams );

	void
	do_set_xyz( core::pose::Pose const & pose, core::Size const i, core::pose::Pose & scratch_pose, core::Size const i_scratch, core::kinematics::Stub const & stub );

	void
	generate_beta_database_test();

	utility::vector1< std::string > load_s_and_l();

} //protein
} //rotamer_sampler
} //protocols

#endif
