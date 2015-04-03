// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/JumpStepWiseSampler.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_JumpStepWiseSampler_HH
#define INCLUDED_protocols_sampler_JumpStepWiseSampler_HH

#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>
#include <protocols/stepwise/sampler/JumpStepWiseSampler.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/Jump.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace stepwise {
namespace sampler {

	class JumpStepWiseSampler: public StepWiseSamplerSized {

	public:

		//constructor
		JumpStepWiseSampler( Size const which_jump,
								 utility::vector1< core::kinematics::Jump > const & jumps,
								 bool const choose_random = false );

		//constructor
		JumpStepWiseSampler();

		//destructor
		~JumpStepWiseSampler();

	public:

		/// @brief Get the total number of rotamers in sampler
		virtual core::Size size() const { return jumps_.size(); }

		/// @brief Apply the i-th rotamer to pose
		virtual void apply( core::pose::Pose &, core::Size const );

		/// @brief Name of the class
		virtual std::string get_name() const { return "JumpStepWiseSampler"; }

		/// @brief Type of class (see enum in StepWiseSamplerTypes.hh)
		virtual StepWiseSamplerType type() const { return JUMP; }

	protected:
		Size which_jump_;
		utility::vector1< core::kinematics::Jump > jumps_;

	};

} //sampler
} //stepwise
} //protocols

#endif
