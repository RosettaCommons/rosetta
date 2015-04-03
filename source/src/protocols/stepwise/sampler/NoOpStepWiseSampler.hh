// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/NoOpStepWiseSampler.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_NoOpStepWiseSampler_HH
#define INCLUDED_protocols_sampler_NoOpStepWiseSampler_HH

#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>
#include <protocols/stepwise/sampler/NoOpStepWiseSampler.fwd.hh>

namespace protocols {
namespace stepwise {
namespace sampler {

	class NoOpStepWiseSampler: public sampler::StepWiseSamplerSized {

	public:

		//constructor
		NoOpStepWiseSampler();

		//destructor
		~NoOpStepWiseSampler();

	public:

		/// @brief Get the total number of rotamers in sampler
		virtual core::Size size() const{ return 1;}

		/// @brief Apply the i-th rotamer to pose
		virtual void apply( core::pose::Pose&, core::Size const ){}

		/// @brief Name of the class
		virtual std::string get_name() const { return "NoOpStepWiseSampler"; }

		/// @brief Type of class (see enum in StepWiseSamplerTypes.hh)
		virtual StepWiseSamplerType type() const { return NO_OP; }

	};

} //sampler
} //stepwise
} //protocols

#endif
