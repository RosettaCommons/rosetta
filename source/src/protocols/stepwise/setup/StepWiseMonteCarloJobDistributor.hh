// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/setup/StepWiseMonteCarloJobDistributor.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_setup_StepWiseMonteCarloJobDistributor_HH
#define INCLUDED_protocols_stepwise_setup_StepWiseMonteCarloJobDistributor_HH

#include <protocols/stepwise/setup/StepWiseJobDistributor.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.fwd.hh>
#include <protocols/stepwise/setup/StepWiseMonteCarloJobDistributor.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace setup {

class StepWiseMonteCarloJobDistributor: public StepWiseJobDistributor {

public:

	//constructor
	StepWiseMonteCarloJobDistributor( stepwise::monte_carlo::StepWiseMonteCarloOP stepwise_monte_carlo,
		std::string const & silent_file,
		core::Size const nstruct );

	//destructor
	~StepWiseMonteCarloJobDistributor();

public:

	virtual std::string get_name() const {
		return "StepWiseMonteCarloJobDistributor";
	}

	virtual
	void
	apply( core::pose::Pose & pose );

	virtual
	void
	initialize( core::pose::Pose const & pose );

	virtual
	bool
	has_another_job();

private:

	void
	move_forward_to_next_model();

	bool
	get_out_tag();

private:

	core::Size count_;
	std::string out_tag_;

	bool init_tag_is_done_;
	std::map< std::string, bool > tag_is_done_;


};

} //setup
} //stepwise
} //protocols

#endif
