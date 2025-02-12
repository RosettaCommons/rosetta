// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/jump/JumpStepWiseSampler.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_JumpStepWiseSampler_HH
#define INCLUDED_protocols_sampler_JumpStepWiseSampler_HH

#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>
#include <protocols/stepwise/sampler/jump/JumpStepWiseSampler.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/Jump.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace jump {

class JumpStepWiseSampler: public StepWiseSamplerSized {

public:

	//constructor
	JumpStepWiseSampler( core::Size const which_jump,
		utility::vector1< core::kinematics::Jump > const & jumps,
		bool const choose_random = false );

	//constructor
	JumpStepWiseSampler();

	//destructor
	~JumpStepWiseSampler() override;

public:

	/// @brief Get the total number of rotamers in sampler
	core::Size size() const override { return jumps_.size(); }

	/// @brief Apply the i-th rotamer to pose
	void apply( core::pose::Pose &, core::Size const ) override;

	/// @brief Name of the class
	std::string get_name() const override { return "JumpStepWiseSampler"; }

	/// @brief Type of class (see enum in toolbox::SamplerPlusPlusTypes.hh)
	toolbox::SamplerPlusPlusType type() const override { return toolbox::JUMP; }

protected:
	core::Size which_jump_;
	utility::vector1< core::kinematics::Jump > jumps_;

};

} //jump
} //sampler
} //stepwise
} //protocols

#endif
