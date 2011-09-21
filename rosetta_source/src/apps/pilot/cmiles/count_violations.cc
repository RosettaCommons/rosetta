// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/cmiles/count_violations.cc
/// @author Christopher Miles (cmiles@uw.edu)

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>

// Project headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>

const std::string ERROR_PREFIX = "Failed to specify required option ";

static basic::Tracer TR("count_violations");

void check_required_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  // structure
  if (!option[in::file::s].user())
    utility_exit_with_message(ERROR_PREFIX + "in:file:s");

  // constraints
  if (!option[constraints::cst_file].user())
    utility_exit_with_message(ERROR_PREFIX + "constraints:cst_file");
}

int main(int argc, char* argv[]) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::Size;
  using core::pose::Pose;
  using core::scoring::ScoreFunction;
  using core::scoring::ScoreFunctionFactory;
  using core::scoring::constraints::ConstraintSetCOP;

  devel::init(argc, argv);
  check_required_options();

  // structure
  Pose pose = *core::import_pose::pose_from_pdb(option[in::file::s]()[1]);

  // score function
  ScoreFunction score = *ScoreFunctionFactory::create_score_function("score0");
  score.set_weight(core::scoring::vdw, 0);

  // constraints
  core::scoring::constraints::add_constraints_from_cmdline(pose, score);
  ConstraintSetCOP constraint_set = pose.constraint_set();

  Size num_violations = constraint_set->show_violations(TR.Trace, pose, 100);
  TR << "Number of constraints violated: " << num_violations << std::endl;
}
