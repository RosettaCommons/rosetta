// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/StepWiseSampleAndScreen.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_StepWiseSampleAndScreen_HH
#define INCLUDED_protocols_stepwise_StepWiseSampleAndScreen_HH

#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/StepWiseSampleAndScreen.fwd.hh>
#include <protocols/stepwise/screener/StepWiseScreener.fwd.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerBase.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace stepwise {

class StepWiseSampleAndScreen: public utility::pointer::ReferenceCount {

public:

	//constructor
	StepWiseSampleAndScreen( sampler::StepWiseSamplerBaseOP sampler,
		utility::vector1< screener::StepWiseScreenerOP > screener );

	//destructor
	~StepWiseSampleAndScreen() override;


public:

	void
	run();

	void
	reset();

	void
	output_counts() const;

	Size const &
	num_tries() const;

	Size const &
	num_successes() const;

	Size
	num_screeners() const;

	void
	output_info_on_random_trials() const;

	void set_max_ntries( core::Size const & setting ){ max_ntries_ = setting; }
	core::Size max_ntries() const{ return max_ntries_; }

	void set_num_random_samples( core::Size const & setting ){ num_random_samples_ = setting; }
	core::Size num_random_samples() const{ return num_random_samples_; }

	void set_verbose( bool const & setting ){ verbose_ = setting; }
	bool verbose() const{ return max_ntries_; }

private:

	void set_ok_to_increment();

	void early_exit_check( Size const n );

public:

	sampler::StepWiseSamplerBaseOP sampler_;
	utility::vector1< screener::StepWiseScreenerOP > screeners_;
	core::Size max_ntries_;
	core::Size num_random_samples_;
	bool verbose_;

};

} //stepwise
} //protocols

#endif
