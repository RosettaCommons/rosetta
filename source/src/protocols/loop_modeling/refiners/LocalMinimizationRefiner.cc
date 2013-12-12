// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/refiners/LocalMinimizationRefiner.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>

// Utility headers
#include <utility/vector1.hh>
#include <iostream>
#include <ctime>

namespace protocols {
namespace loop_modeling {
namespace refiners {

using namespace std;
using core::pose::Pose;
using core::scoring::ScoreFunctionCOP;
using protocols::loops::Loop;

LocalMinimizationRefiner::LocalMinimizationRefiner() {
	minimizer_ = NULL;
}

void LocalMinimizationRefiner::setup(
		Pose & pose, Loop const & loop, ScoreFunctionOP score_function) {

	using core::pose::symmetry::is_symmetric;
	using protocols::simple_moves::MinMover;
	using protocols::simple_moves::symmetry::SymMinMover;
	using protocols::loops::loop_mover::loops_set_chainbreak_weight;

	// The original code checks to see if a certain (non-differentiable) term is 
	// present, and explicitly ignores it if it is.  I took this out because it 
	// required copying the score function, which would make ramping the score 
	// function weights more of a challenge.

	minimizer_ = is_symmetric(pose) ? new SymMinMover() : new MinMover();
	minimizer_->min_type("dfpmin");
	minimizer_->tolerance(1e-3);
	minimizer_->nb_list(true);
	minimizer_->deriv_check(false);
	minimizer_->cartesian(false);

	loops_set_chainbreak_weight(score_function, 1);
	
	apply(pose, loop, score_function);
}

bool LocalMinimizationRefiner::apply(
		Pose & pose, Loop const & loop, ScoreFunctionCOP score_function) {

	using core::kinematics::MoveMapOP;
	using protocols::loops::move_map_from_loops;

	MoveMapOP move_map = move_map_from_loop(pose, loop, false, 10.0, false);

	minimizer_->score_function(score_function);
	minimizer_->movemap(move_map);
	minimizer_->apply(pose);

	return true;
}

}
}
}


