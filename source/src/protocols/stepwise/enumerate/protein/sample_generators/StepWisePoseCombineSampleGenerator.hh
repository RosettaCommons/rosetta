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


#ifndef INCLUDED_protocols_stepwise_StepWisePoseCombineSampleGenerator_HH
#define INCLUDED_protocols_stepwise_StepWisePoseCombineSampleGenerator_HH

#include <protocols/stepwise/enumerate/protein/sample_generators/StepWisePoseSampleGenerator.hh>
#include <protocols/stepwise/enumerate/protein/InputStreamWithResidueInfo.fwd.hh>
#include <protocols/stepwise/enumerate/protein/StepWisePoseSetup.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {
namespace sample_generators {

	class StepWisePoseCombineSampleGenerator: public StepWisePoseSampleGenerator {
	public:

    StepWisePoseCombineSampleGenerator( utility::vector1< InputStreamWithResidueInfoOP > const & input_streams	);

    StepWisePoseCombineSampleGenerator( utility::vector1< InputStreamWithResidueInfoOP > const & input_streams	, StepWisePoseSetupOP stepwise_pose_setup );

		void reset();

		bool has_another_sample();

		void get_next_sample( core::pose::Pose & pose );

		Size size() const;

		void
		set_stepwise_pose_setup( StepWisePoseSetupOP stepwise_pose_setup );

	private:

		utility::vector1< InputStreamWithResidueInfoOP > input_streams_;
		bool need_to_initialize_;

		StepWisePoseSetupOP stepwise_pose_setup_;

  };

} //sample_generators
} //protein
} //enumerate
} //stepwise
} //protocols

#endif

