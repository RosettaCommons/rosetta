// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/BrokenBase.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/BrokenBase.hh>

// C/C++ headers
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/nonlocal.OptionKeys.gen.hh>
#include <numeric/util.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <protocols/abinitio/MaxSeqSepConstraintSet.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/SaneMinMover.hh>

namespace protocols {
namespace nonlocal {

static basic::Tracer TR("protocols.nonlocal.BrokenBase");

void BrokenBase::apply(core::pose::Pose& pose) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::scoring::constraints::ConstraintSetOP;
  using protocols::abinitio::MaxSeqSepConstraintSet;
  using protocols::abinitio::MaxSeqSepConstraintSetOP;

  // Add applicable cutpoint variants and constraints
  add_cutpoint_variants(&pose);
  core::scoring::constraints::add_constraints_from_cmdline_to_pose(pose);

  // Take note of the current constraints, as we are responsible for restoring them
  ConstraintSetOP orig_constraints = pose.constraint_set()->clone();

  // Initialize a MaxSeqSepConstraintSet with the current set of constraints
  MaxSeqSepConstraintSetOP new_constraints = new MaxSeqSepConstraintSet(*orig_constraints, pose.fold_tree());
  new_constraints->set_max_seq_sep(pose.total_residue());
  pose.constraint_set(new_constraints);

  // Apply each of the movers in order
  for (core::Size i = 1; i <= stage_movers_.size(); ++i) {
    show_stage_header(i, TR);
    stage_movers_[i]->apply(pose);

    if (option[OptionKeys::abinitio::debug]())
      emit_intermediate(pose, i);
  }

  // Remove cutpoint variants and restore the pose's original constraints
  remove_cutpoint_variants(&pose);
  pose.constraint_set(orig_constraints);
}

void BrokenBase::add_mover(MoverOP mover) {
  stage_movers_.push_back(mover);
}

void BrokenBase::emit_intermediate(const core::pose::Pose& pose, core::Size stage_num) const {
  std::stringstream ss;
  ss << "nla_stage_" << stage_num << ".pdb";
  core::io::pdb::dump_pdb(pose, ss.str());
}

void BrokenBase::add_cutpoint_variants(core::pose::Pose* pose) const {
  const core::kinematics::FoldTree& tree(pose->fold_tree());
  for (core::Size i = 1; i <= pose->total_residue(); ++i) {
    if (!tree.is_cutpoint(i) || i >= (pose->total_residue() - 1))
      continue;

    core::pose::add_variant_type_to_pose_residue(*pose, core::chemical::CUTPOINT_LOWER, i);
    core::pose::add_variant_type_to_pose_residue(*pose, core::chemical::CUTPOINT_UPPER, i+1);
  }
}

void BrokenBase::remove_cutpoint_variants(core::pose::Pose* pose) const {
  const core::kinematics::FoldTree& tree(pose->fold_tree());
  for (core::Size i = 1; i <= pose->total_residue(); ++i) {
    if (!tree.is_cutpoint(i) || i >= (pose->total_residue() - 1))
      continue;

    core::pose::remove_variant_type_from_pose_residue(*pose, core::chemical::CUTPOINT_LOWER, i);
    core::pose::remove_variant_type_from_pose_residue(*pose, core::chemical::CUTPOINT_UPPER, i+1);
  }
}

void BrokenBase::show_stage_header(core::Size stage_num, std::ostream& out) const {
  out << std::endl << "== Stage " << stage_num << " ==" << std::endl;
}

/// @detail An arbitrary, but hopefully intuitive, coarse-grained ramping strategy.
core::scoring::ScoreFunctionOP BrokenBase::score_function(int stage, int num_residues) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::scoring::ScoreFunctionOP;
  using core::scoring::ScoreFunctionFactory;
  using core::scoring::methods::EnergyMethodOptions;

  TR << "Stage " << stage << " settings:" << std::endl;
  ScoreFunctionOP score = ScoreFunctionFactory::create_score_function("score3");

  // linear_chainbreak ramped up: y=log2(x).
  double m = option[OptionKeys::jumps::increase_chainbreak]();
  score->set_weight(core::scoring::linear_chainbreak, m * (2 + numeric::log(stage, 2)));
  TR << "linear_chainbreak => "
     << score->get_weight(core::scoring::linear_chainbreak)
     << std::endl;

  // atom_pair_constraint and coordinate_constraint ramped down: y=log10(x) / 2.
  if (option[constraints::cst_file].user()) {
    double n = option[OptionKeys::constraints::increase_constraints]();
    score->set_weight(core::scoring::atom_pair_constraint, n * (option[OptionKeys::constraints::cst_weight]() - std::log10((double)stage) / 2));
    TR << "atom_pair_constraint => "
       << score->get_weight(core::scoring::atom_pair_constraint) << std::endl;

    score->set_weight(core::scoring::coordinate_constraint, n * option[OptionKeys::constraints::cst_weight]());
    TR << "coordinate_constraint => "
       << score->get_weight(core::scoring::coordinate_constraint) << std::endl;
  }

  // Residual dipolar coupling held constant
  if (option[OptionKeys::in::file::rdc].user()) {
    score->set_weight(core::scoring::rdc, option[OptionKeys::nonlocal::rdc_weight]());
    TR << "rdc => "
       << score->get_weight(core::scoring::rdc)
       << std::endl;
  }

  // Allowable sequence separation ramped up: y=1.5^x.
  EnergyMethodOptions options(score->energy_method_options());
  double separation = (8 + std::pow(1.5, stage)) * num_residues / 100.0;
  options.cst_max_seq_sep(static_cast<core::Size>(separation));
  score->set_energy_method_options(options);
  TR << "maximum sequence separation => "
     << separation
     << std::endl
     << std::endl;

  // Perturb score function weights with low probability
  if (numeric::random::uniform() < option[OptionKeys::abinitio::prob_perturb_weights]())
    score->perturb_weights();

  return score;
}

// TODO(cmiles) consider how to define a more optimizable gradient
protocols::moves::MoverOP BrokenBase::make_minimizer(core::scoring::ScoreFunctionOP score) {
  using core::kinematics::MoveMap;
  using core::kinematics::MoveMapOP;
  using core::optimization::MinimizerOptions;
  using core::optimization::MinimizerOptionsOP;
  using protocols::moves::SaneMinMover;
  assert(score);

  // define minimizable degrees of freedom
  MoveMapOP movable = new MoveMap();
  movable->set_bb(false);
  movable->set_chi(false);
  movable->set_jump(true);

  // minimizer options
  MinimizerOptionsOP options = new MinimizerOptions("dfpmin", 1e-20, false, false, false);
  options->nblist_auto_update(true);
  options->max_iter(10e5);

  return new SaneMinMover(movable, score, options);
}

}  // namespace nonlocal
}  // namespace protocols
