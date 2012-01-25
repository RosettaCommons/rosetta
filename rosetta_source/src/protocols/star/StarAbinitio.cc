// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/star/StarAbinitio.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/star/StarAbinitio.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <numeric/interpolate.hh>
#include <numeric/prob_util.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/id/SequenceMapping.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/util/ChainbreakUtil.hh>
#include <core/util/kinematics_util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/constraints_additional/MaxSeqSepConstraintSet.hh>
#include <protocols/comparative_modeling/ThreadingJob.hh>
#include <protocols/comparative_modeling/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>
#include <protocols/medal/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/nonlocal/BiasedFragmentMover.hh>
#include <protocols/nonlocal/Policy.hh>
#include <protocols/nonlocal/PolicyFactory.hh>
#include <protocols/nonlocal/SingleFragmentMover.hh>
#include <protocols/nonlocal/util.hh>
#include <protocols/simple_moves/rational_mc/RationalMonteCarlo.hh>

// Package headers
#include <protocols/star/Extender.hh>

namespace protocols {
namespace star {

static basic::Tracer TR("protocols.star.StarAbinitio");

typedef utility::vector1<double> Probabilities;

void compute_per_residue_probabilities(unsigned num_residues,
                                       unsigned fragment_len,
                                       const protocols::loops::Loops& aligned,
                                       const core::kinematics::FoldTree& tree,
                                       Probabilities* probs) {
  assert(probs);
  probs->resize(num_residues, 1);

  // Sampling probability proportional to distance from nearest cut
  Probabilities p_cut;
  protocols::medal::cutpoint_probabilities(num_residues, tree, &p_cut);

  // Lower sampling probability near termini
  Probabilities p_end;
  protocols::medal::end_bias_probabilities(num_residues, &p_end);

  // Product of probabilities
  numeric::product(probs->begin(), probs->end(), p_cut.begin(), p_cut.end());
  numeric::product(probs->begin(), probs->end(), p_end.begin(), p_end.end());

  // Zero-out probabilities of aligned residues and preceding torsions that
  // would cause them to move
  for (unsigned i = 1; i <= aligned.num_loop(); ++i) {
    unsigned stop = aligned[i].stop();
    unsigned start = aligned[i].start() - fragment_len + 1;
    if (start < 1) {
      start = 1;
    }

    for (unsigned j = start; j <= stop; ++j) {
      (*probs)[j] = 0;
    }
  }
  numeric::normalize(probs->begin(), probs->end());

  // Zero-out probabilities of residues that would allow folding across the cut
  protocols::medal::invalidate_residues_spanning_cuts(tree, fragment_len, probs);
  numeric::normalize(probs->begin(), probs->end());
  numeric::print_probabilities(*probs, TR.Debug);
}

void update_chainbreak(core::scoring::ScoreFunctionOP score, double weight) {
  assert(score);
  assert(weight >= 0);
  score->set_weight(core::scoring::linear_chainbreak, weight);
}

/// @detail Utility method for creating a score function of the given type,
/// adding constraints, and updating the linear_chainbreak term
core::scoring::ScoreFunctionOP setup_score(const std::string& weights, double cb) {
  using core::scoring::ScoreFunctionOP;
  using core::scoring::ScoreFunctionFactory;

  ScoreFunctionOP score = ScoreFunctionFactory::create_score_function(weights);
  core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn(*score);
  update_chainbreak(score, cb);
  return score;
}

/// @detail Regulates the application of constraints during folding based on
/// distance between residues in the fold tree. The MonteCarlo object should
/// be reset after calling this function.
void update_sequence_separation(unsigned distance, core::pose::Pose* pose) {
  assert(pose);
  using protocols::constraints_additional::MaxSeqSepConstraintSet;
  using protocols::constraints_additional::MaxSeqSepConstraintSetOP;

  // Regulate application of constraints based on distance between residues in the fold tree
  MaxSeqSepConstraintSetOP new_cst = new MaxSeqSepConstraintSet(*pose->constraint_set(), pose->fold_tree());
  new_cst->set_max_seq_sep(distance);
  pose->constraint_set(new_cst);
  TR << "max_seq_sep => " << distance << std::endl;
}

/// @detail Creates a sequence separation dependent constraint set from restraints
/// specified on the command line and adds it to pose. Takes no action if restraints
/// were not provided. Initial sequence separation threshold is 0.
void setup_constraints(core::pose::Pose* pose) {
  assert(pose);

  // Reads constraints from command line and adds them to pose
  core::scoring::constraints::add_constraints_from_cmdline_to_pose(*pose);
  update_sequence_separation(0, pose);
}

/// @detail If provided, removes restraints from pose.
void tear_down_constraints(core::pose::Pose* pose) {
  assert(pose);
  pose->remove_constraints();
}

/// @detail Writes pose to disk if debug mode enabled
void emit_intermediate(const core::pose::Pose& pose, const std::string& filename) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  if (option[OptionKeys::abinitio::debug]()) {
    pose.dump_pdb(filename);
  }
}

void to_centroid(core::pose::Pose* pose) {
  assert(pose);
  if (!pose->is_centroid()) {
    core::util::switch_to_residue_type_set(*pose, core::chemical::CENTROID);
  }
}

/// @detail Restores simple kinematics to pose
void simple_fold_tree(core::pose::Pose* pose) {
  using core::kinematics::FoldTree;
  assert(pose);
  pose->fold_tree(FoldTree(pose->total_residue()));
}

/// @detail Configures the specified RationalMonteCarlo instance
void configure(protocols::moves::MoverOP mover,
               core::scoring::ScoreFunctionOP score,
               unsigned cycles,
               double temp,
               bool recover,
               protocols::simple_moves::rational_mc::RationalMonteCarlo* rmc) {
  using core::pose::Pose;
  rmc->set_mover(mover);
  rmc->set_temperature(temp);
  rmc->set_num_trials(cycles);
  rmc->set_recover_low(recover);

  // Replacing the score function triggers a reevaluation of cached low/last poses,
  // which fails if none exist.
  const Pose& low = rmc->lowest_score_pose();
  const Pose& last = rmc->last_accepted_pose();

  if (low.total_residue() && last.total_residue()) {
    rmc->set_score_function(score);
  }
}

void StarAbinitio::setup_kinematics(const protocols::loops::Loops& aligned,
                                    const utility::vector1<unsigned>& interior_cuts,
                                    core::pose::Pose* pose) const {
  using core::kinematics::FoldTree;
  using numeric::xyzVector;
  using protocols::loops::Loop;
  using protocols::loops::Loops;
  using utility::vector1;

  assert(pose);
  assert(aligned.num_loop() >= 2);
  assert(interior_cuts.size() == (aligned.num_loop() - 1));

  const unsigned num_residues = pose->total_residue();
  const unsigned vres = num_residues + 1;

  xyzVector<double> center;
  aligned.center_of_mass(*pose, &center);
  core::pose::addVirtualResAsRoot(center, *pose);

  vector1<std::pair<int, int> > jumps;
  for (unsigned i = 1; i <= aligned.num_loop(); ++i) {
    jumps.push_back(std::make_pair(vres, aligned[i].midpoint()));
  }

  vector1<int> cuts(interior_cuts);
  cuts.push_back(num_residues);

  ObjexxFCL::FArray2D_int ft_jumps(2, jumps.size());
  for (unsigned i = 1; i <= jumps.size(); ++i) {
    ft_jumps(1, i) = std::min(jumps[i].first, jumps[i].second);
    ft_jumps(2, i) = std::max(jumps[i].first, jumps[i].second);
  }

  ObjexxFCL::FArray1D_int ft_cuts(cuts.size());
  for (unsigned i = 1; i <= cuts.size(); ++i) {
    ft_cuts(i) = cuts[i];
  }

  FoldTree tree(vres);
  bool status = tree.tree_from_jumps_and_cuts(vres,          // nres_in
                                              jumps.size(),  // num_jump_in
                                              ft_jumps,      // jump_point
                                              ft_cuts,       // cuts
                                              vres);         // root

  assert(status);
  pose->fold_tree(tree);
  core::util::add_cutpoint_variants(pose);

  TR << pose->fold_tree() << std::endl;
}

void StarAbinitio::apply(core::pose::Pose& pose) {
  using core::kinematics::ShortestPathInFoldTree;
  using core::scoring::Energies;
  using core::scoring::ScoreFunctionOP;
  using protocols::comparative_modeling::ThreadingJob;
  using protocols::loops::Loops;
  using protocols::loops::LoopsOP;
  using protocols::moves::MoverOP;
  using protocols::nonlocal::BiasedFragmentMover;
  using protocols::nonlocal::PolicyOP;
  using protocols::nonlocal::PolicyFactory;
  using protocols::simple_moves::rational_mc::RationalMonteCarlo;
  using utility::vector1;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  ThreadingJob const * const job = protocols::nonlocal::current_job();
  assert(job);

  to_centroid(&pose);
  emit_intermediate(pose, "star_initial.pdb");

  Extender extender(job->alignment().clone(), pose.total_residue());
  extender.set_secondary_structure(pred_ss_);
  extender.extend_unaligned(&pose);
  emit_intermediate(pose, "star_extended.pdb");

  const Loops& aligned = *(extender.aligned());
  const Loops& unaligned = *(extender.unaligned());
  TR << "Aligned: " << aligned << std::endl;
  TR << "Unaligned: " << unaligned << std::endl;

  // Prepare the pose
  const unsigned num_residues = pose.total_residue();
  setup_kinematics(aligned, extender.cutpoints(), &pose);
  setup_constraints(&pose);

  Probabilities probs_sm, probs_lg;
  compute_per_residue_probabilities(num_residues, fragments_sm_->max_frag_length(), aligned, pose.fold_tree(), &probs_sm);
  compute_per_residue_probabilities(num_residues, fragments_lg_->max_frag_length(), aligned, pose.fold_tree(), &probs_lg);

  // Simulation parameters
  const int num_stages = 4;
  const double t = option[OptionKeys::abinitio::temperature]();
  const double m = option[OptionKeys::abinitio::increase_cycles]();
  const double c = option[OptionKeys::jumps::increase_chainbreak]();

  // Compute ramping schedules
  ShortestPathInFoldTree ft_dist(pose.fold_tree());
  int max_dist = ft_dist.max_dist();

  vector1<int> seq_sep(num_stages);
  for (int i = 1; i <= num_stages; ++i) {
    seq_sep[i] = static_cast<int>(numeric::linear_interpolate<int>(0, max_dist, i, num_stages));
  }

  // Movers
  MoverOP fragments_lg_uni = new BiasedFragmentMover(PolicyFactory::get_policy("uniform", fragments_lg_), probs_lg);
  MoverOP fragments_sm_uni = new BiasedFragmentMover(PolicyFactory::get_policy("uniform", fragments_sm_), probs_sm);
  MoverOP fragments_sm_smo = new BiasedFragmentMover(PolicyFactory::get_policy("smooth", fragments_sm_), probs_sm);

  // Scores
  ScoreFunctionOP score1  = setup_score("score0", c * 0.10);
  ScoreFunctionOP score2  = setup_score("score1", c * 0.25);
  ScoreFunctionOP score3a = setup_score("score2", c * 0.50);
  ScoreFunctionOP score3b = setup_score("score5", c * 0.75);
  ScoreFunctionOP score4  = setup_score("score3", c * 1.00);

  // Stage 1
  TR << "Stage 1" << std::endl;
  update_sequence_separation(seq_sep[1], &pose);

  RationalMonteCarlo stage_mover(fragments_lg_uni, score1, static_cast<unsigned>(m * 2000), t, false);  // 1a
  stage_mover.apply(pose);

  configure(fragments_sm_uni, score1, static_cast<unsigned>(m * 2000), t, false, &stage_mover);  // 1b
  stage_mover.apply(pose);

  emit_intermediate(pose, "star_stage_1.pdb");

  // Stage 2
  TR << "Stage 2" << std::endl;
  update_sequence_separation(seq_sep[2], &pose);

  configure(fragments_lg_uni, score2, static_cast<unsigned>(m * 4000), t, true, &stage_mover);  // 2a
  stage_mover.apply(pose);

  configure(fragments_sm_uni, score2, static_cast<unsigned>(m * 4000), t, true, &stage_mover);  // 2b
  stage_mover.apply(pose);

  emit_intermediate(pose, "star_stage_2.pdb");

  // Stage 3
  TR << "Stage 3" << std::endl;
  update_sequence_separation(seq_sep[3], &pose);
  for (unsigned i = 1; i <= 10; ++i) {
    ScoreFunctionOP score = ((i % 2) == 0 && i <= 7) ? score3a : score3b;

    configure(fragments_lg_uni, score, static_cast<unsigned>(m * 4000), t, true, &stage_mover);  // 3a
    stage_mover.apply(pose);

    configure(fragments_sm_uni, score, static_cast<unsigned>(m * 4000), t, true, &stage_mover);  // 3b
    stage_mover.apply(pose);
  }
  emit_intermediate(pose, "star_stage_3.pdb");

  // Stage 4
  TR << "Stage 4" << std::endl;
  update_sequence_separation(seq_sep[4], &pose);
  configure(fragments_sm_smo, score4, static_cast<unsigned>(m * 8000), t, true, &stage_mover);
  for (unsigned i = 1; i <= 3; ++i) {
    stage_mover.apply(pose);  // 4
  }
  emit_intermediate(pose, "star_stage_4.pdb");

  // Removing the virtual residue introduced by the star fold tree invalidates
  // the pose's cached energies. Doing so causes the score line in the silent
  // file to be 0.
  Energies energies = pose.energies();

  // Housekeeping
  tear_down_constraints(&pose);
  tear_down_kinematics(&pose);
  pose.set_new_energies_object(new Energies(energies));
}

void StarAbinitio::tear_down_kinematics(core::pose::Pose* pose) const {
  assert(pose);

  core::Size vres = pose->total_residue();
  pose->conformation().delete_residue_slow(vres);

  simple_fold_tree(pose);
  core::util::remove_cutpoint_variants(pose);
}

StarAbinitio::StarAbinitio() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::fragment::FragmentIO;
  using core::fragment::FragSetOP;
  using core::fragment::SecondaryStructure;

  FragmentIO io;
  fragments_lg_ = io.read_data(option[in::file::frag9]());
  fragments_sm_ = io.read_data(option[in::file::frag3]());

  // Approximate secondary structure from fragments when psipred isn't available
  if (option[in::file::psipred_ss2].user()) {
    pred_ss_ = new SecondaryStructure();
    pred_ss_->read_psipred_ss2(option[in::file::psipred_ss2]());
  } else {
    pred_ss_ = new SecondaryStructure(*fragments_sm_);
  }
}

std::string StarAbinitio::get_name() const {
  return "StarAbinitio";
}

protocols::moves::MoverOP StarAbinitio::clone() const {
  return new StarAbinitio(*this);
}

protocols::moves::MoverOP StarAbinitio::fresh_instance() const {
  return new StarAbinitio();
}

}  // namespace star
}  // namespace protocols
