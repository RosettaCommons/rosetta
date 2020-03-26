// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerSized.hh
/// @brief Abstract Base Class for Sampler sampler with finite size.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_stepwise_sampler_StepWiseSamplerSized_HH
#define INCLUDED_stepwise_sampler_StepWiseSamplerSized_HH

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerSized.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSampler.hh>

namespace protocols {
namespace stepwise {
namespace sampler {

class StepWiseSamplerSized: public sampler::StepWiseSampler {

public:
	StepWiseSamplerSized();

	~StepWiseSamplerSized() override;

	/// @brief Get the total number of rotamers in sampler
	virtual core::Size size() const = 0;

	/// @brief Initialization
	void init() override;

	/// @brief Reset to the first (or random if random()) rotamer.
	void reset() override;

	/// @brief Move to next rotamer
	void operator++() override;

	/// @brief Check if reach the end of rotamer list
	bool not_end() const override;

	/// @brief Apply the i-th rotamer to pose
	virtual void apply( core::pose::Pose&, core::Size const ) = 0;

	/// @brief Apply the current rotamer to pose
	void apply( core::pose::Pose & pose ) override { apply( pose, id_ ); }

	/// @brief Name of the class
	std::string get_name() const override { return "StepWiseSamplerSized"; }

	/// @brief Type of class (see enum in SamplerPlusPlusTypes.hh)
	toolbox::SamplerPlusPlusType type() const override { return toolbox::SIZED; }

	core::Size const & id() const { return id_; }

	virtual void fast_forward(){ id_ = size(); }

	virtual void set_id( core::Size const setting ){ id_ = setting; }

	/// @brief output summary of class
	void show( std::ostream & out, core::Size const indent = 0 ) const override {
		for ( core::Size n = 1; n <= indent; n++ ) out << ' ';
		out << get_name() << " [" << size() << ']' << std::endl;
	}

protected:

	core::Size id_;

};

} //sampler
} //stepwise
} //protocols
#endif
