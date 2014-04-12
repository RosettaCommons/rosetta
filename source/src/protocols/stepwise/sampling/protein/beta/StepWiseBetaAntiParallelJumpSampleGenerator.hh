// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ResidueSampler.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das
/// @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_stepwise_protein_StepWiseBetaAntiParallelJumpSampleGenerator_HH
#define INCLUDED_protocols_stepwise_protein_StepWiseBetaAntiParallelJumpSampleGenerator_HH

#include <protocols/stepwise/sampling/protein/sample_generator/StepWiseProteinJumpSampleGenerator.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <core/kinematics/Jump.hh>

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

	class StepWiseBetaAntiParallelJumpSampleGenerator: public StepWiseProteinJumpSampleGenerator {
	public:

    StepWiseBetaAntiParallelJumpSampleGenerator( core::pose::Pose const & pose, core::Size const moving_residue );


		core::Size
		get_antiparallel_beta_jumps( core::pose::Pose const & pose,
																 int const sample_res );

	private:
		utility::vector1< core::kinematics::Jump > jumps_;

  };

} //protein
} //sampling
} //stepwise
} //protocols

#endif

