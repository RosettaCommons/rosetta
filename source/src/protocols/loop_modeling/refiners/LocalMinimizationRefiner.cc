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
#include <core/kinematics/MoveMap.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
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

bool LocalMinimizationRefiner::do_apply(Pose & pose) {
	using core::kinematics::MoveMapOP;
	using core::pose::symmetry::is_symmetric;
	using protocols::loops::move_map_from_loops;
	using protocols::simple_moves::MinMover;
	using protocols::simple_moves::MinMoverOP;
	using protocols::simple_moves::symmetry::SymMinMover;

	pose.update_residue_neighbors();

	Loops const & loops = get_loops();
	ScoreFunctionCOP score_function = get_score_function();
	MoveMapOP move_map = move_map_from_loops(pose, loops, false, 10.0, false);
	MinMoverOP minimizer = is_symmetric(pose) ? new SymMinMover : new MinMover;

	minimizer->min_type("dfpmin");
	minimizer->tolerance(1e-3);
	minimizer->nb_list(true);
	minimizer->deriv_check(false);
	minimizer->cartesian(false);
	minimizer->score_function(score_function);
	minimizer->movemap(move_map);
	minimizer->apply(pose);

	return true;
}

}
}
}


