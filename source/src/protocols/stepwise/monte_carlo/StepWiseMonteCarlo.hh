// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/monte_carlo/StepWiseMonteCarlo.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_StepWiseMonteCarlo_HH
#define INCLUDED_protocols_stepwise_monte_carlo_StepWiseMonteCarlo_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.fwd.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>
#include <protocols/stepwise/monte_carlo/submotif/SubMotifLibrary.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <core/pose/Pose.fwd.hh>
#include <ctime>

namespace protocols {
namespace stepwise {
namespace monte_carlo {

class StepWiseMonteCarlo: public protocols::moves::Mover {

public:

	//constructor
	StepWiseMonteCarlo( core::scoring::ScoreFunctionCOP scorefxn );

	//destructor
	~StepWiseMonteCarlo();

public:

	virtual std::string get_name() const {
		return "StepWiseMonteCarlo";
	}

	/// @brief Apply the loop-rebuild protocol to the input pose
	virtual
	void apply ( core::pose::Pose & pose );

	/// @brief setter for native poses contained for rms ---- we should get rid of this method? it is widely used, but a bit unsafe
	virtual void set_native_pose( core::pose::PoseCOP pose );

	void
	set_options( options::StepWiseMonteCarloOptionsCOP options );

	options::StepWiseMonteCarloOptionsCOP
	options() const { return options_; }

	void set_model_tag( std::string const & setting ){ model_tag_ = setting; }
	std::string model_tag() const{ return model_tag_; }

	void set_out_path( std::string const & setting ){ out_path_ = setting; }
	std::string out_path() const{ return out_path_; }

	void set_move( mover::StepWiseMove const & setting ){ move_ = setting; }
	mover::StepWiseMove move() const { return move_; }

	void set_submotif_library( monte_carlo::submotif::SubMotifLibraryCOP setting );

	mover::StepWiseMasterMoverCOP master_mover() const { return master_mover_; }

private:

	void initialize();

	void initialize_scorefunction();

	void initialize_for_movie( core::pose::Pose const & pose );

	void do_main_loop( core::pose::Pose & pose );

	void
	output_movie( core::pose::Pose const & pose, core::Size const k, std::string const & tag, std::string const & movie_file );

	core::Real
	display_progress( core::pose::Pose & pose, core::Size const cycle_num );

	core::Real show_scores( core::pose::Pose & pose, std::string const & tag );

	void
	anneal_missing( protocols::moves::MonteCarloOP monte_carlo );

private:

	core::scoring::ScoreFunctionCOP scorefxn_input_;

	core::scoring::ScoreFunctionOP scorefxn_;
	options::StepWiseMonteCarloOptionsCOP options_;
	mover::StepWiseMasterMoverOP master_mover_;

	// are these really in use anymore?
	core::Real const max_missing_weight_;
	core::Real missing_weight_interval_;
	core::Real missing_weight_;

	// for movies
	std::string model_tag_;
	std::string out_path_;
	std::string movie_file_trial_;
	std::string movie_file_accepted_;

	// for testing individual moves
	mover::StepWiseMove move_;

	// timing poses
	std::clock_t start_time_;
};

} //monte_carlo
} //stepwise
} //protocols

#endif
