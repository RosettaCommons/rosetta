// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/medal/MedalExchangeMover.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/medal/MedalExchangeMover.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>

// Project headers
#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/util/SwitchResidueTypeSet.hh>
// AUTO-REMOVED #include <protocols/abinitio/MaxSeqSepConstraintSet.hh>
// AUTO-REMOVED #include <protocols/comparative_modeling/util.hh>
#include <protocols/jd2/ThreadingJob.hh>
// AUTO-REMOVED #include <protocols/loops/LoopRelaxMover.hh>
#include <protocols/loops/LoopRelaxThreadingMover.hh>
#include <protocols/loops/Loops.hh>
// AUTO-REMOVED #include <protocols/loops/util.hh>
#include <protocols/nonlocal/StarTreeBuilder.hh>
#include <protocols/nonlocal/util.hh>

// Package headers
// AUTO-REMOVED #include <protocols/medal/util.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace medal {

static basic::Tracer TR("protocols.medal.MedalExchangeMover");

MedalExchangeMover::MedalExchangeMover() : MedalMover() {}

void MedalExchangeMover::apply(core::pose::Pose& pose) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::scoring::ScoreFunctionOP;
  using core::scoring::methods::EnergyMethodOptions;
  using protocols::jd2::ThreadingJob;
  using protocols::loops::LoopRelaxThreadingMover;
  using protocols::loops::Loops;
  using namespace protocols::moves;
  using namespace protocols::nonlocal;

  ThreadingJob const * const job = protocols::nonlocal::current_job();
  const unsigned num_residues = pose.total_residue();

  // Build up a threading model
  LoopRelaxThreadingMover closure;
  closure.setup();
  closure.apply(pose);

  // Configure the score functions used in the simulation
  core::util::switch_to_residue_type_set(pose, core::chemical::CENTROID);
  core::scoring::constraints::add_constraints_from_cmdline_to_pose(pose);

  // Decompose the structure into chunks
  Loops chunks;
  decompose_structure(pose, &chunks);

  StarTreeBuilder builder;
  builder.set_up(chunks, &pose);
  TR << pose.fold_tree() << std::endl;

  // Compute per-residue sampling probabilities
  Probabilities probs;
  compute_per_residue_probabilities(num_residues, job->alignment(), chunks, pose.fold_tree(), *fragments_sm_, &probs);

  // Fragment insertion with short range constraints and no chainbreak
  protocols::nonlocal::add_cutpoint_variants(&pose);
  ScoreFunctionOP score = score_function();
  EnergyMethodOptions options(score->energy_method_options());

  options.cst_max_seq_sep(option[OptionKeys::rigid::short_range_seqsep]());
  score->set_energy_method_options(options);
  MoverOP stage1 = create_fragment_mover(score, fragments_sm_, probs, "uniform", 25);
  stage1->apply(pose);

  // Fragment insertion with medium range constraints and no chainbreak
  options.cst_max_seq_sep(option[OptionKeys::rigid::medium_range_seqsep]());
  score->set_energy_method_options(options);
  MoverOP stage2 = create_fragment_mover(score, fragments_sm_, probs, "uniform", 25);
  stage2->apply(pose);

  // Fragment insertion with long range constraints and chainbreak
  options.cst_max_seq_sep(num_residues);
  score->set_energy_method_options(options);
  score->set_weight(core::scoring::linear_chainbreak, 2 * option[OptionKeys::jumps::increase_chainbreak]());
  MoverOP stage3 = create_fragment_and_rigid_mover(pose, score, fragments_sm_, probs, "uniform", 25);
  stage3->apply(pose);

  // Loop closure
  builder.tear_down(&pose);
  do_loop_closure(&pose);
  protocols::nonlocal::remove_cutpoint_variants(&pose);

  // Return to centroid representation and rescore
  core::util::switch_to_residue_type_set(pose, core::chemical::CENTROID);
  score_pose(*score, "Final model", &pose);
}

std::string MedalExchangeMover::get_name() const {
  return "MedalExchangeMover";
}

protocols::moves::MoverOP MedalExchangeMover::clone() const {
  return new MedalExchangeMover(*this);
}

protocols::moves::MoverOP MedalExchangeMover::fresh_instance() const {
  return new MedalExchangeMover();
}

}  // namespace medal
}  // namespace protocols
