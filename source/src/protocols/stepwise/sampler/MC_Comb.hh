// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/MC_Comb.hh
/// @brief Ensemble of Markov chain samplers for modeler combinatorially.
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_sampler_MC_Comb_HH
#define INCLUDED_protocols_sampler_MC_Comb_HH

// Unit headers
#include <protocols/stepwise/sampler/MC_Comb.fwd.hh>

// Package headers
#include <protocols/stepwise/sampler/MC_StepWiseSampler.hh>

namespace protocols {
namespace stepwise {
namespace sampler {

class MC_Comb : public MC_StepWiseSampler {
public:
	MC_Comb();

	virtual ~MC_Comb(){}

	/// @brief Initialization
	virtual void init();

	/// @brief Reset to the first (or random if random()) rotamer
	virtual void reset();

	/// @brief Move to next rotamer
	virtual void operator++();

	/// @brief Update the DOFs
	virtual void update();

	/// @brief Apply the current rotamer to pose
	virtual void apply( core::pose::Pose & pose );

	/// @brief Get uniform modeler state
	virtual void set_uniform_modeler( bool const setting );

	/// @brief Add one more rotamer sampler to this sampler
	virtual void add_external_loop_rotamer( MC_StepWiseSamplerOP const & rotamer ) {
		rotamer_list_.push_back( rotamer );
		set_init( false );
	}

	/// @brief Clear all rotamer samplers stored in this sampler
	virtual void clear_rotamer() {
		rotamer_list_.clear();
		set_init( false );
	}

	/// @brief Name of the class
	virtual std::string get_name() const { return "MC_Comb"; }

private:
	bool is_empty_;
	utility::vector1<MC_StepWiseSamplerOP> rotamer_list_;
};

} //sampler
} //stepwise
} //protocols

#endif

