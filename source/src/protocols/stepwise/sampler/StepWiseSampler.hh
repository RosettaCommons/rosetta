// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/StepWiseSampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampler_StepWiseSampler_HH
#define INCLUDED_protocols_stepwise_sampler_StepWiseSampler_HH

#include <protocols/toolbox/SamplerPlusPlus.hh>
#include <protocols/stepwise/sampler/StepWiseSampler.fwd.hh>

namespace protocols {
namespace stepwise {
namespace sampler {

class StepWiseSampler: public toolbox::SamplerPlusPlus {

public:

	//constructor
	StepWiseSampler():
		SamplerPlusPlus(),
		random_( true )
	{}

	//destructor
	~StepWiseSampler()
	{}

public:

	/// @brief Check if there are more rotamers available
	virtual bool not_end() const { return true; }

	/// @brief Check if is random modeler
	virtual bool random() const { return random_; }

	/// @brief Set the random modeler state
	virtual void set_random( bool const setting ) { random_ = setting; }

private:

	bool random_;

};

} //sampler
} //stepwise
} //protocols

#endif
