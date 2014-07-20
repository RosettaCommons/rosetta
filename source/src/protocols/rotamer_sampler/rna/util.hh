// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/rna/util.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_rna_RNA_RotamerSamplerSamplerUtil_HH
#define INCLUDED_protocols_rotamer_sampler_rna_RNA_RotamerSamplerSamplerUtil_HH

#include <protocols/rotamer_sampler/RotamerSamplerBase.fwd.hh>
#include <protocols/stepwise/sampling/modeler_options/StepWiseModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

using namespace core;

namespace protocols {
namespace rotamer_sampler {
namespace rna {

	rotamer_sampler::RotamerSamplerBaseOP
	setup_rotamer_sampler( core::pose::Pose const & pose,
												 protocols::stepwise::sampling::modeler_options::StepWiseModelerOptionsCOP options,
												 protocols::stepwise::sampling::working_parameters::StepWiseWorkingParametersCOP working_parameters,
												 bool const build_pose_from_scratch,
												 bool const kic_sampling,
												 bool const close_chain );

	bool
	sampling_sugar_at_five_prime( pose::Pose const & pose,
																Size const moving_suite );


	bool
	sampling_sugar_at_three_prime( pose::Pose const & pose,
																 Size const moving_suite );


} //rna
} //rotamer_sampler
} //protocols

#endif
