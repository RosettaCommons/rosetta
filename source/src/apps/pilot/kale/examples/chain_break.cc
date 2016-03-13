// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Protocols headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>

// Utility headers
#include <devel/init.hh>
#include <utility/vector1.hh>
#include <iostream>
#include <ctime>

using namespace std;
using core::pose::Pose;
using core::import_pose::pose_from_file;
using core::kinematics::FoldTree;
using core::scoring::ScoreFunctionOP;
using protocols::loops::Loop;
using protocols::loops::Loops;
using protocols::loops::fold_tree_from_loops;

int main(int argc, char** argv) {
	devel::init(argc, argv);

	Pose pose; pose_from_file(pose, "chain_break_input.pdb", core::import_pose::PDB_file);
	Loop loop(36, 47, 47);
	Loops loops; loops.add_loop(loop);
	ScoreFunctionOP score_function = core::scoring::get_score_function();
	FoldTree tree; protocols::loops::fold_tree_from_loops(pose, loops, tree);
	pose.fold_tree(tree);

	using core::pose::symmetry::is_symmetric;
  using protocols::simple_moves::MinMover;
  using protocols::simple_moves::MinMoverOP;
  using protocols::loops::loop_mover::loops_set_chainbreak_weight;

  MinMoverOP minimizer = new MinMover();
  minimizer->min_type("lbfgs_armijo_nonmonotone");
  minimizer->tolerance(1e-3);
  minimizer->nb_list(true);
  minimizer->deriv_check(false);
  minimizer->cartesian(false);

  loops_set_chainbreak_weight(score_function, 1);

  using core::kinematics::MoveMapOP;
  using protocols::loops::move_map_from_loops;

  MoveMapOP move_map = move_map_from_loop(pose, loop, false, 10.0, false);

  minimizer->score_function(score_function);
  minimizer->movemap(move_map);
  minimizer->apply(pose);

  pose.dump_pdb("chain_break_output.pdb");
}
