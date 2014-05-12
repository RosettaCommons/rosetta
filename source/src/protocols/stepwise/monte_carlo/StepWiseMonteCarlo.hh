// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/StepWiseMonteCarlo.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_StepWiseMonteCarlo_HH
#define INCLUDED_protocols_stepwise_monte_carlo_StepWiseMonteCarlo_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarlo.fwd.hh>
#include <protocols/stepwise/monte_carlo/StepWiseMonteCarloOptions.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/AddMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/DeleteMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/FromScratchMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/AddOrDeleteMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/ResampleMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/SWA_Move.hh>
#include <protocols/stepwise/sampling/StepWiseModeler.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <core/pose/Pose.hh>

/*
using namespace core;

Commented out because “using namespace X” in header files outside of class declaration is explicitly forbidden
by our coding convention due to problems it create on modern compilers and because of the name clashing.
For more information please see: https://wiki.rosettacommons.org/index.php/Coding_conventions#Using
*/


namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace rna {

	class StepWiseMonteCarlo: public protocols::moves::Mover {

	public:

	//constructor
		StepWiseMonteCarlo( core::scoring::ScoreFunctionCOP scorefxn );

		//destructor
		~StepWiseMonteCarlo();

		virtual std::string get_name() const {
			return "StepWiseMonteCarlo";
		}

		/// @brief Apply the loop-rebuild protocol to the input pose
		virtual
		void apply ( core::pose::Pose & pose );

	public:

		void
		set_options( StepWiseMonteCarloOptionsCOP options );

		void set_model_tag( std::string const & setting ){ model_tag_ = setting; }
		std::string model_tag() const{ return model_tag_; }

		void set_out_path( std::string const & setting ){ out_path_ = setting; }
		std::string out_path() const{ return out_path_; }

		void set_move( SWA_Move const setting ){ move_ = setting; }

		void set_enumerate( bool const setting ){ enumerate_ = setting; }

		void set_just_preminimize( bool const setting ){	just_preminimize_ = setting; }

	private:

		void initialize_scorefunction();

		void initialize_movers();

		void initialize_for_movie( pose::Pose const & pose );

		void do_main_loop( pose::Pose & pose );

		void initialize_pose_if_empty( pose::Pose & pose );

		void
		output_movie( pose::Pose const & pose, Size const k, std::string const tag, std::string const & movie_file );

		Real
		display_progress( pose::Pose & pose, Size const cycle_num );

		std::string
		get_all_res_list( pose::Pose & pose );

		Real show_scores( core::pose::Pose & pose, std::string const tag );

		void
		set_minimize_single_res( bool const minimize_single_res );

		void
		anneal_missing( protocols::moves::MonteCarloOP monte_carlo );

		bool
		do_test_move( pose::Pose & pose );

		sampling::StepWiseModelerOP
		setup_unified_stepwise_modeler();

		void
		preminimize_pose( pose::Pose & pose );

	private:

		core::scoring::ScoreFunctionCOP scorefxn_input_;
		core::scoring::ScoreFunctionOP scorefxn_;

		StepWiseMonteCarloOptionsCOP options_;
		sampling::StepWiseModelerOP stepwise_modeler_;
		DeleteMoverOP delete_mover_;
		AddMoverOP add_mover_;
		FromScratchMoverOP from_scratch_mover_;
		AddOrDeleteMoverOP add_or_delete_mover_;
		ResampleMoverOP resample_mover_;

		bool minimize_single_res_;
		core::Real const max_missing_weight_;
		core::Real missing_weight_interval_;
		core::Real missing_weight_;

		// for movies
		std::string model_tag_;
		std::string out_path_;
		std::string movie_file_trial_;
		std::string movie_file_accepted_;

		// for testing individual moves
		SWA_Move move_;
		bool enumerate_;
		bool just_preminimize_;

	};

} //rna
} //monte_carlo
} //stepwise
} //protocols

#endif
