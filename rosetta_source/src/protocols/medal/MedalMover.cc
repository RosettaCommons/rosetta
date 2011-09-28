// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/medal/MedalMover.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/medal/MedalMover.hh>

// C/C++ headers
#include <iostream>
#include <string>

// External headers
#include <boost/unordered/unordered_map.hpp>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/abinitio/MaxSeqSepConstraintSet.hh>
#include <protocols/comparative_modeling/util.hh>
#include <protocols/jd2/ThreadingJob.hh>
#include <protocols/loops/LoopRelaxMover.hh>
#include <protocols/loops/LoopRelaxThreadingMover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>
#include <protocols/moves/RationalMonteCarlo.hh>
#include <protocols/moves/RigidBodyMotionMover.hh>
#include <protocols/nonlocal/Policy.hh>
#include <protocols/nonlocal/PolicyFactory.hh>
#include <protocols/nonlocal/SingleFragmentMover.hh>
#include <protocols/nonlocal/StarTreeBuilder.hh>
#include <protocols/nonlocal/util.hh>

namespace protocols {
namespace medal {

typedef boost::unordered_map<int, core::kinematics::Jump> Jumps;

static basic::Tracer TR("protocols.medal.MedalMover");

void MedalMover::apply(core::pose::Pose& pose) {
  using core::scoring::ScoreFunctionOP;
  using protocols::jd2::ThreadingJob;
  using protocols::loops::LoopRelaxThreadingMover;
  using protocols::loops::Loops;
  using protocols::nonlocal::StarTreeBuilder;

  ThreadingJob const * const job = protocols::nonlocal::current_job();

  // Build up a threading model
  LoopRelaxThreadingMover closure;
  closure.setup();
  closure.apply(pose);

  // Decompose the structure into chunks based on consecutive CA-CA distances
  Loops chunks;
  protocols::nonlocal::chunks_by_CA_CA_distance(pose, &chunks);

  // Configure the score functions used in the simulation
  core::util::switch_to_residue_type_set(pose, core::chemical::CENTROID);
  ScoreFunctionOP perturb_score = perturb_score_function();
  ScoreFunctionOP insertion_score = insert_score_function();

  // Add constraints to the pose, score the initial model
  core::scoring::constraints::add_constraints_from_cmdline_to_pose(pose);
  perturb_score->show(TR, pose);

  // Define the kinematics
  StarTreeBuilder builder;
  builder.set_up(chunks, &pose);
  TR << pose.fold_tree() << std::endl;

  // Action:
  //  1. Rigid body perturbation
  //  2. Fragment insertion (uniform policy)
  //  3. CCD closure
  protocols::nonlocal::add_cutpoint_variants(&pose);
  do_rigid_body_moves(perturb_score, &pose);
  do_fragment_insertion(insertion_score, &pose);

  // Housekeeping
  builder.tear_down(&pose);
  protocols::nonlocal::remove_cutpoint_variants(&pose);

  // Close remaining chainbreaks
  do_loop_closure(&pose);
}

void MedalMover::do_loop_closure(core::pose::Pose* pose) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::fragment::FragmentIO;
  using core::fragment::FragSetOP;
  using core::kinematics::FoldTree;
  using protocols::loops::LoopRelaxMover;
  using protocols::loops::Loops;
  assert(pose);

  // Choose chainbreaks automatically
  Loops empty;
  LoopRelaxMover closure;
  closure.remodel("quick_ccd");
  closure.intermedrelax("no");
  closure.refine("no");
  closure.relax("no");
  closure.loops(empty);

  // Small fragments
  FragmentIO io;
  FragSetOP fragments_sm = io.read_data(option[in::file::frag3]());
  utility::vector1<core::fragment::FragSetOP> fragments;
  fragments.push_back(fragments_sm);
  closure.frag_libs(fragments);

  // Use atom pair constraints when available
  closure.cmd_line_csts(option[constraints::cst_fa_file].user());

  // Simple kinematics
  FoldTree tree(pose->total_residue());
  pose->fold_tree(tree);
  closure.apply(*pose);

  core::util::switch_to_residue_type_set(*pose, core::chemical::CENTROID);
}

void MedalMover::do_rigid_body_moves(const core::scoring::ScoreFunctionOP& score,
                                     core::pose::Pose* pose) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using protocols::moves::MoverOP;
  using protocols::moves::RationalMonteCarlo;
  using protocols::moves::RigidBodyMotionMover;
  assert(pose);

  Jumps jumps;
  jumps_from_pose(*pose, &jumps);

  MoverOP rigid_mover = new RationalMonteCarlo(
      new RigidBodyMotionMover(jumps),
      score,
      option[OptionKeys::rigid::rigid_body_cycles](),
      option[OptionKeys::rigid::temperature](),
      true);

  rigid_mover->apply(*pose);
}

void MedalMover::do_fragment_insertion(const core::scoring::ScoreFunctionOP& score,
                                       core::pose::Pose* pose) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::Size;
  using core::fragment::FragmentIO;
  using core::fragment::FragSetOP;
  using core::kinematics::MoveMap;
  using core::kinematics::MoveMapOP;
  using protocols::moves::MoverOP;
  using protocols::moves::RationalMonteCarlo;
  using protocols::nonlocal::PolicyOP;
  using protocols::nonlocal::PolicyFactory;
  using protocols::nonlocal::SingleFragmentMover;
  assert(pose);

  FragmentIO io;
  FragSetOP fragments = io.read_data(option[in::file::frag3]());

  // Fragment selection policies
  const Size num_fragments = option[OptionKeys::abinitio::number_3mer_frags]();
  PolicyOP uniform = PolicyFactory::get_policy("uniform", fragments, num_fragments);
  PolicyOP smooth  = PolicyFactory::get_policy("smooth", fragments, num_fragments);

  // Define degrees of freedom
  MoveMapOP movable = new MoveMap();
  movable->set_bb(true);

  // Small fragments, uniform selection
  MoverOP uniform_mover = new RationalMonteCarlo(
      new SingleFragmentMover(fragments, movable, uniform),
      score,
      option[OptionKeys::rigid::fragment_cycles](),
      option[OptionKeys::rigid::temperature](),
      true);

  uniform_mover->apply(*pose);

  // Small fragments, smooth selection
  MoverOP smooth_mover = new RationalMonteCarlo(
      new SingleFragmentMover(fragments, movable, smooth),
      score,
      option[OptionKeys::rigid::fragment_cycles](),
      option[OptionKeys::rigid::temperature](),
      true);

  smooth_mover->apply(*pose);
}

void MedalMover::jumps_from_pose(const core::pose::Pose& pose, Jumps* jumps) const {
  using core::kinematics::Jump;
  assert(jumps);

  for (core::Size i = 1; i <= pose.num_jump(); ++i) {
    const Jump& jump = pose.jump(i);
    (*jumps)[i] = jump;
    TR.Debug << "Added jump_num " << i << ": " << jump << std::endl;
  }
}

core::scoring::ScoreFunctionOP MedalMover::base_score_function() const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::scoring::ScoreFunctionFactory;
  using core::scoring::ScoreFunctionOP;
  using core::scoring::methods::EnergyMethodOptions;

  ScoreFunctionOP score = ScoreFunctionFactory::create_score_function("score0");

  // Maximum sequence separation
  EnergyMethodOptions options(score->energy_method_options());
  options.cst_max_seq_sep(option[OptionKeys::rigid::sequence_separation]());
  score->set_energy_method_options(options);

  // Enable specific energy terms
  score->set_weight(core::scoring::atom_pair_constraint, option[OptionKeys::constraints::cst_weight]());
  score->set_weight(core::scoring::hbond_lr_bb, 1);
  score->set_weight(core::scoring::hbond_sr_bb, 1);
  score->set_weight(core::scoring::linear_chainbreak, option[OptionKeys::jumps::increase_chainbreak]());
  score->set_weight(core::scoring::rg, option[OptionKeys::abinitio::rg_reweight]());
  score->set_weight(core::scoring::sheet, 1);
  score->set_weight(core::scoring::vdw, 0.25);

  // Disable specific energy terms
  score->set_weight(core::scoring::rama, 0);

  core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score);
  return score;
}

// TODO(cmiles) specialize
core::scoring::ScoreFunctionOP MedalMover::perturb_score_function() const {
  return base_score_function();
}

// TODO(cmiles) specialize
core::scoring::ScoreFunctionOP MedalMover::insert_score_function() const {
  return base_score_function();
}

std::string MedalMover::get_name() const {
  return "MedalMover";
}

protocols::moves::MoverOP MedalMover::clone() const {
  return new MedalMover(*this);
}

protocols::moves::MoverOP MedalMover::fresh_instance() const {
  return new MedalMover();
}

}  // namespace medal
}  // namespace protocols
