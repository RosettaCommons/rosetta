// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/BrokenFold.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/BrokenFold.hh>

// C/C++ headers
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

// External headers
#include <boost/format.hpp>

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
#include <protocols/moves/RationalMonteCarlo.hh>
#include <protocols/moves/SaneMinMover.hh>

// Package headers
#include <protocols/nonlocal/Policy.hh>
#include <protocols/nonlocal/PolicyFactory.hh>
#include <protocols/nonlocal/SingleFragmentMover.hh>
#include <protocols/nonlocal/util.hh>

namespace protocols {
namespace nonlocal {

static basic::Tracer TR("protocols.nonlocal.BrokenFold");

BrokenFold::BrokenFold(core::fragment::FragSetOP fragments_lg,
                       core::fragment::FragSetOP fragments_sm,
                       core::kinematics::MoveMapOP movable) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::Real;
  using core::Size;
  using protocols::moves::MoverOP;
  using protocols::moves::RationalMonteCarlo;

  // Consider only the top k best fragments in each library
  Size num_fragments = option[OptionKeys::abinitio::number_3mer_frags];
  Size nres = fragments_sm->max_pos();

  // Fragment selection policies
  PolicyOP uniform_lg = PolicyFactory::get_policy("uniform", fragments_lg, num_fragments);
  PolicyOP uniform_sm = PolicyFactory::get_policy("uniform", fragments_sm, num_fragments);
  PolicyOP smooth_sm  = PolicyFactory::get_policy("smooth", fragments_sm, num_fragments);

  // Bookkeeping
  Size stage = 1;
  Real temperature = option[OptionKeys::abinitio::temperature]();

  // Large fragments, uniformly chosen
  for (Size i = 1; i <= 3; ++i, ++stage) {
    MoverOP minimize = minimizer(score_function(stage, nres));
    MoverOP fold = new RationalMonteCarlo(
        new SingleFragmentMover(fragments_lg, movable, uniform_lg),
        score_function(stage, nres),
        cycles(stage),
        temperature,
        true);

    movers_.push_back(minimize);
    movers_.push_back(fold);
  }

  // Small fragments, uniformly chosen
  for (Size i = 1; i <= 3; ++i, ++stage) {
    MoverOP minimize = minimizer(score_function(stage, nres));
    MoverOP fold = new RationalMonteCarlo(
        new SingleFragmentMover(fragments_sm, movable, uniform_sm),
        score_function(stage, nres),
        cycles(stage),
        temperature,
        true);

    movers_.push_back(minimize);
    movers_.push_back(fold);
  }

  // Small fragments, smoothly chosen
  for (Size i = 1; i <= 3; ++i, ++stage) {
    MoverOP minimize = minimizer(score_function(stage, nres));
    MoverOP fold = new RationalMonteCarlo(
        new SingleFragmentMover(fragments_sm, movable, smooth_sm),
        score_function(stage, nres),
        cycles(stage),
        temperature,
        true);

    movers_.push_back(minimize);
    movers_.push_back(fold);
  }
}

void BrokenFold::add_cutpoint_variants(core::pose::Pose* pose) const {
  const core::kinematics::FoldTree& tree(pose->fold_tree());
  for (core::Size i = 1; i <= pose->total_residue(); ++i) {
    if (!tree.is_cutpoint(i) || i >= (pose->total_residue() - 1))
      continue;

    core::pose::add_variant_type_to_pose_residue(*pose, core::chemical::CUTPOINT_LOWER, i);
    core::pose::add_variant_type_to_pose_residue(*pose, core::chemical::CUTPOINT_UPPER, i+1);
  }
}

void BrokenFold::remove_cutpoint_variants(core::pose::Pose* pose) const {
  const core::kinematics::FoldTree& tree(pose->fold_tree());
  for (core::Size i = 1; i <= pose->total_residue(); ++i) {
    if (!tree.is_cutpoint(i) || i >= (pose->total_residue() - 1))
      continue;

    core::pose::remove_variant_type_from_pose_residue(*pose, core::chemical::CUTPOINT_LOWER, i);
    core::pose::remove_variant_type_from_pose_residue(*pose, core::chemical::CUTPOINT_UPPER, i+1);
  }
}

void BrokenFold::show_stage_header(core::Size stage_num, std::ostream& out) const {
  out << std::endl << "== Stage " << stage_num << " ==" << std::endl;
}

/// @detail An arbitrary, but hopefully intuitive, coarse-grained ramping strategy.
core::scoring::ScoreFunctionOP BrokenFold::score_function(int stage, int num_residues) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::scoring::ScoreFunctionOP;
  using core::scoring::ScoreFunctionFactory;
  using core::scoring::methods::EnergyMethodOptions;

  ScoreFunctionOP score = ScoreFunctionFactory::create_score_function("score3");

  // linear_chainbreak ramped up: y=log2(x).
  double m = option[OptionKeys::jumps::increase_chainbreak]();
  score->set_weight(core::scoring::linear_chainbreak, m * (2 + numeric::log(stage, 2)));

  // atom_pair_constraint ramped down: y=log10(x) / 2.
  // coordinate_constraint held constant.
  if (option[constraints::cst_file].user()) {
    double n = option[OptionKeys::constraints::increase_constraints]();
    score->set_weight(core::scoring::atom_pair_constraint, n * (option[OptionKeys::constraints::cst_weight]() - std::log10((double)stage) / 2));
    score->set_weight(core::scoring::coordinate_constraint, n * option[OptionKeys::constraints::cst_weight]());
  }

  // Residual dipolar coupling held constant
  if (option[OptionKeys::in::file::rdc].user())
    score->set_weight(core::scoring::rdc, option[OptionKeys::nonlocal::rdc_weight]());

  // Allowable sequence separation ramped up: y=1.5^x.
  EnergyMethodOptions options(score->energy_method_options());
  double separation = (8 + std::pow(1.5, stage)) * num_residues / 100.0;
  options.cst_max_seq_sep(static_cast<core::Size>(separation));
  score->set_energy_method_options(options);

  // Perturb score function weights with low probability
  if (numeric::random::uniform() < option[OptionKeys::abinitio::prob_perturb_weights]())
    score->perturb_weights();

  return score;
}

void BrokenFold::apply(core::pose::Pose& pose) {
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
  for (core::Size i = 1; i <= movers_.size(); ++i) {
    show_stage_header(i, TR);
    movers_[i]->apply(pose);

    if (option[OptionKeys::abinitio::debug]())
      emit_intermediate(pose, str(boost::format("nla_stage_%d.pdb") % i));
  }

  // Remove cutpoint variants and restore the pose's original constraints
  remove_cutpoint_variants(&pose);
  pose.constraint_set(orig_constraints);
}

protocols::moves::MoverOP BrokenFold::minimizer(core::scoring::ScoreFunctionOP score) {
  using core::kinematics::MoveMap;
  using core::kinematics::MoveMapOP;
  using core::optimization::MinimizerOptions;
  using core::optimization::MinimizerOptionsOP;
  using protocols::moves::SaneMinMover;
  assert(score);

  score->set_weight(core::scoring::chainbreak, 1);
  score->set_weight(core::scoring::linear_chainbreak, 1);
  score->set_weight(core::scoring::overlap_chainbreak, 1);
  score->set_weight(core::scoring::distance_chainbreak, 1);

  // define minimizable degrees of freedom
  MoveMapOP movable = new MoveMap();
  movable->set_bb(false);
  movable->set_chi(false);
  movable->set_jump(true);

  // minimizer options
  MinimizerOptionsOP options = new MinimizerOptions("dfpmin", 1e-10, true);
  options->nblist_auto_update(true);
  options->max_iter(10000);

  return new SaneMinMover(movable, score, options);
}

core::Size BrokenFold::cycles(int stage) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  assert(stage >= 0);

  double m = option[OptionKeys::abinitio::increase_cycles]();
  double cycles = (stage <= 6) ? 2000 * m : 4000 * m;
  return static_cast<core::Size>(cycles);
}

std::string BrokenFold::get_name() const {
  return "BrokenFold";
}

}  // namespace nonlocal
}  // namespace protocols
