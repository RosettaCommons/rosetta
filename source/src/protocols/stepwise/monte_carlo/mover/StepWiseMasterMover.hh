// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_monte_carlo_mover_StepWiseMasterMover_HH
#define INCLUDED_protocols_stepwise_monte_carlo_mover_StepWiseMasterMover_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/AddMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/DeleteMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/FromScratchMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/AddOrDeleteMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/ResampleMover.fwd.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMoveSelector.fwd.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.fwd.hh>
#include <protocols/stepwise/monte_carlo/submotif/SubMotifLibrary.fwd.hh>
#include <protocols/stepwise/modeler/StepWiseModeler.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {

	class StepWiseMasterMover: public protocols::moves::Mover {

	public:

		//constructor
		StepWiseMasterMover( core::scoring::ScoreFunctionCOP scorefxn, options::StepWiseMonteCarloOptionsCOP options );

		//destructor
		~StepWiseMasterMover();

	public:

		using Mover::apply;
		virtual void apply( Pose & );

		virtual std::string get_name() const { return "StepWiseMasterMover"; }

		void
		apply( core::pose::Pose & pose,
					 StepWiseMove const & stepwise_move,
					 bool const figure_out_all_possible_moves = true );

		void
		do_the_move( StepWiseMove const & move, core::pose::Pose & pose );

		bool apply_legacy( Pose & );

		void
		initialize( core::scoring::ScoreFunctionCOP scorefxn, options::StepWiseMonteCarloOptionsCOP options );

		void initialize_pose_if_empty( core::pose::Pose & pose );

		void
		set_minimize_single_res( bool const minimize_single_res );

		bool
		do_test_move( StepWiseMove const & move,
									core::pose::Pose & pose );

		void
		preminimize_pose( core::pose::Pose & pose );

		void
		set_options( options::StepWiseMonteCarloOptionsCOP options ){ options_ = options; }

		void
		set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn ){ scorefxn_ = scorefxn; }

		std::string move_type_string() const { return move_type_string_; }

		bool success() const { return success_; }

		core::Real const & proposal_density_ratio() const { return proposal_density_ratio_; }

		void
		build_full_model( core::pose::Pose const & start_pose, core::pose::Pose & full_model_pose );

		void set_submotif_library( monte_carlo::submotif::SubMotifLibraryCOP setting ) { submotif_library_ = setting; }

	private:

		void initialize();

		bool minimize_single_res_;

		modeler::StepWiseModelerOP
		setup_unified_stepwise_modeler();

		bool
		test_all_moves( core::pose::Pose & pose );

		void
		test_all_moves_recursively( core::pose::Pose & pose );

	private:

		core::scoring::ScoreFunctionCOP scorefxn_;
		options::StepWiseMonteCarloOptionsCOP options_;

		modeler::StepWiseModelerOP stepwise_modeler_;
		DeleteMoverOP delete_mover_;
		AddMoverOP add_mover_;
		FromScratchMoverOP from_scratch_mover_;
		AddOrDeleteMoverOP add_or_delete_mover_;
		ResampleMoverOP resample_mover_;
		monte_carlo::submotif::SubMotifLibraryCOP submotif_library_;
		StepWiseMoveSelectorOP stepwise_move_selector_;

		bool success_;
		std::string move_type_string_;
		core::Real proposal_density_ratio_;

		Size num_tested_moves_;
	};

} //mover
} //monte_carlo
} //stepwise
} //protocols

#endif
