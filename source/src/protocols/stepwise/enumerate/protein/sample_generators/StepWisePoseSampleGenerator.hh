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


#ifndef INCLUDED_protocols_stepwise_StepWisePoseSampleGenerator_HH
#define INCLUDED_protocols_stepwise_StepWisePoseSampleGenerator_HH

#include <core/pose/Pose.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/enumerate/protein/sample_generators/StepWisePoseSampleGenerator.fwd.hh>

//Auto Headers
#include <utility/vector1.hh>


namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {
namespace sample_generators {

	class StepWisePoseSampleGenerator: public utility::pointer::ReferenceCount {
	public:

    StepWisePoseSampleGenerator() {}

    ~StepWisePoseSampleGenerator() {}

		virtual
		void reset() = 0;

		virtual
		bool has_another_sample() = 0;

		virtual
		void get_next_sample( core::pose::Pose & ) = 0;

		virtual
		Size size() const { return 0;}

  };

} //sample_generators
} //protein
} //enumerate
} //stepwise
} //protocols

#endif

