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


#ifndef INCLUDED_protocols_sampler_StepWiseSamplerComb_HH
#define INCLUDED_protocols_sampler_StepWiseSamplerComb_HH

// Unit headers
#include <protocols/stepwise/sampler/StepWiseSamplerComb.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerBase.hh>

namespace protocols {
namespace stepwise {
namespace sampler {

class StepWiseSamplerComb : public StepWiseSamplerBase {
public:
	StepWiseSamplerComb();

	virtual ~StepWiseSamplerComb();

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

	/// @brief Clear all rotamer samplers stored in this sampler
	virtual void clear_rotamer() {
		rotamer_list_.clear();
		set_init( false );
	}

	/// @brief Name of the class
	virtual std::string get_name() const { return "StepWiseSamplerComb"; }

	/// @brief Type of class (see enum in StepWiseSamplerTypes.hh)
	virtual StepWiseSamplerType type() const { return COMB; }

	/// @brief output summary of class
	virtual
	void show( std::ostream & out, Size const indent = 0 ) const;

private:
	bool is_empty_;
	utility::vector1<StepWiseSamplerBaseOP> rotamer_list_;
};

} //sampler
} //stepwise
} //protocols

#endif

