// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/protein/loop_close/util.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_protein_loop_close_util_HH
#define INCLUDED_protocols_stepwise_sampling_protein_loop_close_util_HH

#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/sampling/modeler_options/StepWiseModelerOptions.fwd.hh>
#include <protocols/rotamer_sampler/RotamerSamplerSized.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {
namespace loop_close {

	void
	kic_close_loops_in_samples( rotamer_sampler::RotamerSamplerSizedOP & sampler,
															core::pose::Pose const & pose,
															working_parameters::StepWiseWorkingParametersCOP working_parameters_,
															modeler_options::StepWiseModelerOptionsCOP options_ );

	void
	enable_sampling_of_loop_takeoff( rotamer_sampler::RotamerSamplerSizedOP & sampler,
																	 core::pose::Pose & pose,
																	 working_parameters::StepWiseWorkingParametersCOP working_parameters_,
																	 modeler_options::StepWiseModelerOptionsCOP options_ );


} //loop_close
} //protein
} //sampling
} //stepwise
} //protocols

#endif
