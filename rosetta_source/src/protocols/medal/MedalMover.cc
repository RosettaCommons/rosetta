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
#include <iomanip>
#include <iostream>
#include <string>

// External headers
#include <boost/bind.hpp>
#include <boost/format.hpp>
#include <boost/function.hpp>
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
#include <numeric/prob_util.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/sequence/SequenceAlignment.hh>
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
#include <protocols/nonlocal/BiasedFragmentMover.hh>
#include <protocols/nonlocal/Policy.hh>
#include <protocols/nonlocal/PolicyFactory.hh>
#include <protocols/nonlocal/StarTreeBuilder.hh>
#include <protocols/nonlocal/util.hh>

// Package headers
#include <protocols/medal/util.hh>

namespace protocols {
namespace medal {

typedef boost::function<void(const core::pose::Pose&)> Trigger;
typedef boost::unordered_map<int, core::kinematics::Jump> Jumps;
typedef utility::vector1<double> Probabilities;

static basic::Tracer TR("protocols.medal.MedalMover");

// Assumptions:
//   - There is a 1-to-1 mapping between chunks and jumps
//   - The ordering of chunks is identical locally and in <fragment_mover>
//   - The ordering of jumps is identical locally and in <rigid_body_mover>
void on_pose_accept(const protocols::loops::Loops& chunks,
                    const Jumps& jumps,
                    const core::pose::Pose& pose,
                    protocols::moves::MoverOP rigid_body_mover,
                    protocols::moves::MoverOP fragment_mover) {
  using core::io::silent::SilentFileData;
  using core::io::silent::SilentStructFactory;
  using core::io::silent::SilentStructOP;
  assert(jumps.size() == chunks.size());

  static int num_accepted = 0;

  SilentStructOP silent = SilentStructFactory::get_instance()->get_silent_struct_out();
  silent->fill_struct(pose, str(boost::format("accepted_pose_%d") % num_accepted++));

  SilentFileData sfd;
  sfd.write_silent_struct(*silent, "medal.accepted.out");
}

void compute_per_residue_probabilities(const unsigned num_residues,
                                       const core::sequence::SequenceAlignment& alignment,
                                       const protocols::loops::Loops& chunks,
                                       const core::kinematics::FoldTree& tree,
                                       const core::fragment::FragSet& fragments,
                                       Probabilities* probs) {
  using namespace std;

  Probabilities p_alignment, p_chunk, p_cut, p_end;
  alignment_probabilities(num_residues, alignment, &p_alignment);
  chunk_probabilities(chunks, &p_chunk);
  cutpoint_probabilities(num_residues, tree, &p_cut);
  end_bias_probabilities(num_residues, &p_end);

  // Product of probabilities
  probs->insert(probs->begin(), p_alignment.begin(), p_alignment.end());
  numeric::product(probs->begin(), probs->end(), p_chunk.begin(), p_chunk.end());
  numeric::product(probs->begin(), probs->end(), p_cut.begin(), p_cut.end());
  numeric::product(probs->begin(), probs->end(), p_end.begin(), p_end.end());

  // Zero-out probabilities of residues that would allow folding across the cut
  invalidate_residues_spanning_cuts(tree, fragments.max_frag_length(), probs);
  numeric::normalize(probs->begin(), probs->end());

  // Display individual and combined probabilities
  TR << setw(10) << "Residue"
     << setw(15) << "P(align)"
     << setw(15) << "P(chunk)"
     << setw(15) << "P(cut)"
     << setw(15) << "P(end)"
     << setw(15) << "P(combined)"
     << endl;

  for (unsigned i = 1; i <= probs->size(); ++i) {
    TR << setw(10) << i;
    TR << fixed << setw(15) << setprecision(5) << p_alignment[i];
    TR << fixed << setw(15) << setprecision(5) << p_chunk[i];
    TR << fixed << setw(15) << setprecision(5) << p_cut[i];
    TR << fixed << setw(15) << setprecision(5) << p_end[i];
    TR << fixed << setw(15) << setprecision(5) << (*probs)[i] << endl;
  }
  TR.flush_all_channels();
}

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
  using core::scoring::ScoreFunctionOP;
  using core::sequence::SequenceAlignment;
  using protocols::jd2::ThreadingJob;
  using protocols::loops::LoopRelaxThreadingMover;
  using protocols::loops::Loops;
  using namespace protocols::moves;
  using namespace protocols::nonlocal;

  ThreadingJob const * const job = protocols::nonlocal::current_job();
  const SequenceAlignment& alignment = job->alignment();
  const unsigned num_residues = pose.total_residue();

  // Build up a threading model
  LoopRelaxThreadingMover closure;
  closure.setup();
  closure.apply(pose);

  // Configure the score functions used in the simulation
  core::util::switch_to_residue_type_set(pose, core::chemical::CENTROID);
  core::scoring::constraints::add_constraints_from_cmdline_to_pose(pose);
  ScoreFunctionOP score = score_function();
  score_pose(*score, &pose);

  // Decompose the structure into chunks based on consecutive CA-CA distances
  Loops chunks;
  protocols::nonlocal::chunks_by_CA_CA_distance(pose, &chunks);

  StarTreeBuilder builder;
  builder.set_up(chunks, &pose);
  TR << pose.fold_tree() << std::endl;

  // Compute per-residue sampling probabilities
  Probabilities probs;
  compute_per_residue_probabilities(num_residues, alignment, chunks, pose.fold_tree(), *fragments_sm_, &probs);

  // Rigid body moves
  Jumps jumps;
  core::pose::jumps_from_pose(pose, &jumps);
  MoverOP rigid_body_mover = new RigidBodyMotionMover(jumps);

  // Fragment insertion moves
  MoverOP fragment_mover =
      new BiasedFragmentMover(fragments_sm_,
                              PolicyFactory::get_policy("uniform", fragments_sm_, 25),
                              probs);

  // Alternating rigid body and fragment insertion moves
  CyclicMoverOP meta = new CyclicMover();
  meta->enqueue(rigid_body_mover);
  meta->enqueue(fragment_mover);

  const Size cycles = option[OptionKeys::rigid::rigid_body_cycles]() + option[OptionKeys::rigid::fragment_cycles]();
  RationalMonteCarlo mover(meta, score, cycles, option[OptionKeys::rigid::temperature](), true);

  // Partial function application and callback registration
  Trigger callback = boost::bind(&on_pose_accept, chunks, jumps, _1, rigid_body_mover, fragment_mover);
  mover.add_trigger(callback);

  protocols::nonlocal::add_cutpoint_variants(&pose);
  mover.apply(pose);
  protocols::nonlocal::remove_cutpoint_variants(&pose);

  // Loop closure
  builder.tear_down(&pose);
  do_loop_closure(&pose);

  // Return to centroid representation and rescore
  core::util::switch_to_residue_type_set(pose, core::chemical::CENTROID);
  score_pose(*score, &pose);
}

void MedalMover::score_pose(const core::scoring::ScoreFunction& score,
                            core::pose::Pose* pose) const {
  assert(pose);
  score.show(TR, *pose);
  TR.flush_all_channels();
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

core::scoring::ScoreFunctionOP MedalMover::score_function() const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::scoring::ScoreFunctionFactory;
  using core::scoring::ScoreFunctionOP;
  using core::scoring::methods::EnergyMethodOptions;

  ScoreFunctionOP score = ScoreFunctionFactory::create_score_function("score3");

  // Maximum sequence separation
  EnergyMethodOptions options(score->energy_method_options());
  options.cst_max_seq_sep(option[OptionKeys::rigid::sequence_separation]());
  score->set_energy_method_options(options);

  // Enable specific energy terms
  score->set_weight(core::scoring::atom_pair_constraint, option[OptionKeys::constraints::cst_weight]());
  score->set_weight(core::scoring::hbond_lr_bb, 0.3);
  score->set_weight(core::scoring::hbond_sr_bb, 0.3);
  score->set_weight(core::scoring::linear_chainbreak, option[OptionKeys::jumps::increase_chainbreak]());

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
