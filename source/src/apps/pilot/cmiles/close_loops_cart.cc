// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/cmiles/close_loops_cart.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

int main(int argc, char* argv[]) {
	try {

  using core::kinematics::MoveMap;
  using core::optimization::CartesianMinimizer;
  using core::optimization::MinimizerOptions;
  using core::pose::PoseOP;
  using core::scoring::ScoreFunctionFactory;
  using core::scoring::ScoreFunctionOP;
  using protocols::simple_moves::SwitchResidueTypeSetMover;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  devel::init(argc, argv);

  PoseOP pose = core::import_pose::pose_from_file(option[OptionKeys::in::file::s]()[1], core::import_pose::PDB_file);
  SwitchResidueTypeSetMover to_centroid(core::chemical::CENTROID);
  to_centroid.apply(*pose);
  pose->dump_pdb("starting_cart.pdb");

  ScoreFunctionOP score = ScoreFunctionFactory::create_score_function("score4_smooth_cart");
  core::scoring::constraints::add_constraints_from_cmdline(*pose, *score);
  (*score)(*pose);

  MinimizerOptions options_lbfgs("lbfgs_armijo_nonmonotone", 0.01, true, false, false);
  options_lbfgs.max_iter(1000);

  MoveMap mm;
  mm.set_bb(true);
  mm.set_chi(true);
  mm.set_jump(true);

  CartesianMinimizer minimizer;
  minimizer.run(*pose, mm, *score, options_lbfgs);
  pose->dump_pdb("ending_cart.pdb");

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
