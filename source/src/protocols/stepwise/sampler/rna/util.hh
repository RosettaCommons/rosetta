// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/rna/util.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_rna_RNA_StepWiseSamplerSamplerUtil_HH
#define INCLUDED_protocols_sampler_rna_RNA_StepWiseSamplerSamplerUtil_HH

#include <protocols/stepwise/sampler/StepWiseSamplerBase.fwd.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.fwd.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

sampler::StepWiseSamplerBaseOP
setup_sampler( core::pose::Pose const & pose,
	protocols::stepwise::modeler::options::StepWiseModelerOptionsCOP options,
	protocols::stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP working_parameters,
	bool const build_pose_from_scratch,
	bool const kic_modeler,
	bool const close_chain );

bool
modeler_sugar_at_five_prime( pose::Pose const & pose,
	Size const moving_suite );


bool
modeler_sugar_at_three_prime( pose::Pose const & pose,
	Size const moving_suite );


} //rna
} //sampler
} //stepwise
} //protocols

#endif
