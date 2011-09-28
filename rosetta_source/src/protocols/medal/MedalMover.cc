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
#include <protocols/moves/CyclicMover.hh>
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

MedalMover::MedalMover() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::fragment::FragmentIO;
  using core::fragment::FragSetOP;

  FragmentIO io;
  fragments_lg_ = io.read_data(option[in::file::frag9]());
  fragments_sm_ = io.read_data(option[in::file::frag3]());
}

void MedalMover::apply(core::pose::Pose& pose) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::kinematics::MoveMap;
  using core::kinematics::MoveMapOP;
  using core::scoring::ScoreFunctionOP;
  using protocols::jd2::ThreadingJob;
  using protocols::loops::LoopRelaxThreadingMover;
  using protocols::loops::Loops;
  using protocols::moves::CyclicMover;
  using protocols::moves::CyclicMoverOP;
  using protocols::moves::MoverOP;
  using protocols::moves::RationalMonteCarlo;
  using protocols::moves::RigidBodyMotionMover;
  using protocols::nonlocal::SingleFragmentMover;
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
  ScoreFunctionOP score = score_function();

  // Add constraints to the pose, score the initial model
  core::scoring::constraints::add_constraints_from_cmdline_to_pose(pose);
  score->show(TR, pose);

  // Define the kinematics
  StarTreeBuilder builder;
  builder.set_up(chunks, &pose);
  TR << pose.fold_tree() << std::endl;

  // Define the base movers
  Jumps jumps;
  jumps_from_pose(pose, &jumps);
  MoverOP rigid_body_mover = new RigidBodyMotionMover(jumps);

  MoveMapOP movable = new MoveMap();
  movable->set_bb(true);
  MoverOP fragment_mover = new SingleFragmentMover(fragments_sm_, movable);

  // Alternating rigid body and fragment insertion moves
  CyclicMoverOP meta = new CyclicMover();
  meta->enqueue(rigid_body_mover);
  meta->enqueue(fragment_mover);

  protocols::nonlocal::add_cutpoint_variants(&pose);
  const Size cycles = option[OptionKeys::rigid::rigid_body_cycles]() + option[OptionKeys::rigid::fragment_cycles]();
  RationalMonteCarlo mover(meta, score, cycles, option[OptionKeys::rigid::temperature](), true);
  mover.apply(pose);
  protocols::nonlocal::remove_cutpoint_variants(&pose);

  // Loop closure
  builder.tear_down(&pose);
  do_loop_closure(&pose);

  // Return to centroid representation and rescore
  core::util::switch_to_residue_type_set(pose, core::chemical::CENTROID);
  score->show(TR, pose);
}

void MedalMover::do_loop_closure(core::pose::Pose* pose) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
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

  utility::vector1<core::fragment::FragSetOP> fragments;
  fragments.push_back(fragments_sm_);
  closure.frag_libs(fragments);

  // Use atom pair constraints when available
  closure.cmd_line_csts(option[constraints::cst_fa_file].user());

  // Simple kinematics
  FoldTree tree(pose->total_residue());
  pose->fold_tree(tree);
  closure.apply(*pose);
}

void MedalMover::jumps_from_pose(const core::pose::Pose& pose, Jumps* jumps) const {
  assert(jumps);
  for (core::Size i = 1; i <= pose.num_jump(); ++i) {
    const core::kinematics::Jump& jump = pose.jump(i);
    (*jumps)[i] = jump;
    TR.Debug << "Added jump_num " << i << ": " << jump << std::endl;
  }
}

core::scoring::ScoreFunctionOP MedalMover::score_function() const {
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
  score->set_weight(core::scoring::vdw, 1);

  // Disable specific energy terms
  score->set_weight(core::scoring::rama, 0);

  core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score);
  return score;
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
