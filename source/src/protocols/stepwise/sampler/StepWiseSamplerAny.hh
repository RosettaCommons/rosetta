// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/StepWiseSamplerAny.hh
/// @brief Aggregate multiple samplers for modeler from any one of them.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_sampler_StepWiseSamplerAny_HH
#define INCLUDED_protocols_sampler_StepWiseSamplerAny_HH

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerAny.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerBase.hh>

namespace protocols {
namespace stepwise {
namespace sampler {

class StepWiseSamplerAny : public StepWiseSamplerBase {
public:
	StepWiseSamplerAny();

	virtual ~StepWiseSamplerAny();

	/// @brief Initialization
	virtual void init();

	/// @brief Reset to the first (or random if random()) rotamer
	virtual void reset();

	/// @brief Move to next rotamer
	virtual void operator++();

	/// @brief Check if reach the end of rotamer list
	virtual bool not_end() const;

	/// @brief Apply the current rotamer to pose
	virtual void apply( core::pose::Pose & pose );

	/// @brief Set the random modeler state
	virtual void set_random( bool const setting );

	/// @brief Add one more rotamer sampler to this sampler
	virtual void add_external_loop_rotamer( StepWiseSamplerBaseOP const & rotamer ) {
		rotamer_list_.push_back( rotamer );
		set_init( false );
	}

	/// @brief Add one more rotamer sampler to this sampler
	virtual void add_external_loop_rotamer(
		StepWiseSamplerBaseOP const & rotamer,
		core::Real const weight
	) {
		rotamer_list_.push_back( rotamer );
		weights_.push_back( weight );
		set_init( false );
	}

	/// @brief Set the weights of each rotamer sampler
	virtual void set_weights( utility::vector1<core::Real> const & weights ) {
		weights_ = weights;
		set_init( false );
	}

	/// @brief Clear all rotamer samplers stored in this sampler
	virtual void clear_rotamer() {
		rotamer_list_.clear();
		set_init( false );
	}

	/// @brief Name of the class
	virtual std::string get_name() const { return "StepWiseSamplerAny"; }

	/// @brief Type of class (see enum in StepWiseSamplerTypes.hh)
	virtual StepWiseSamplerType type() const { return ANY; }

private:
	core::Size curr_rotamer_;
	bool is_weighted_, is_empty_, has_empty_;
	utility::vector1<StepWiseSamplerBaseOP> rotamer_list_;
	utility::vector1<core::Real> weights_, cdf_;
};

} //sampler
} //stepwise
} //protocols

#endif

