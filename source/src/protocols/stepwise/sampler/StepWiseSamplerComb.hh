// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerComb.hh
/// @brief Aggregate of multiple rotamer samplers for modeler combinatorially.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_stepwise_sampler_StepWiseSamplerComb_HH
#define INCLUDED_stepwise_sampler_StepWiseSamplerComb_HH

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerComb.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSampler.hh>

namespace protocols {
namespace stepwise {
namespace sampler {

class StepWiseSamplerComb : public sampler::StepWiseSampler {
public:
	StepWiseSamplerComb();

	~StepWiseSamplerComb() override;

	/// @brief Initialization
	void init() override;

	/// @brief Reset to the first (or random if random()) rotamer
	void reset() override;

	/// @brief Move to next rotamer
	void operator++() override;

	/// @brief Check if reach the end of rotamer list
	bool not_end() const override;

	/// @brief Apply the current rotamer to pose
	void apply( core::pose::Pose & pose ) override;

	/// @brief Set the random modeler state
	void set_random( bool const setting ) override;

	/// @brief Add one more rotamer sampler to this sampler
	virtual void add_external_loop_rotamer( sampler::StepWiseSamplerOP const & rotamer ) {
		rotamer_list_.push_back( rotamer );
		set_init( false );
	}

	/// @brief Clear all rotamer samplers stored in this sampler
	virtual void clear_rotamer() {
		rotamer_list_.clear();
		set_init( false );
	}

	/// @brief Name of the class
	std::string get_name() const override { return "StepWiseSamplerComb"; }

	/// @brief Type of class (see enum in SamplerPlusPlusTypes.hh)
	toolbox::SamplerPlusPlusType type() const override { return toolbox::COMB; }

	/// @brief output summary of class
	void show( std::ostream & out, core::Size const indent = 0 ) const override;

private:
	bool is_empty_;
	utility::vector1<sampler::StepWiseSamplerOP> rotamer_list_;
};

} //sampler
} //stepwise
} //protocols

#endif

