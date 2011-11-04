// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/medal/MedalAbinitioMover.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/medal/MedalAbinitioMover.hh>

// C/C++ headers
#include <set>
#include <string>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/nonlocal.OptionKeys.gen.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <numeric/interpolate.hh>
#include <numeric/prob_util.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>
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

typedef std::set<int> Jumps;

MedalAbinitioMover::MedalAbinitioMover() : MedalMover() {}

protocols::moves::MoverOP MedalAbinitioMover::fragment_mover(const core::pose::Pose& pose,
                                                             const core::fragment::FragSetOP fragments,
                                                             const protocols::nonlocal::PolicyOP policy) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::fragment::FragSetOP;
  using protocols::nonlocal::BiasedFragmentMover;
  using utility::vector1;

  // Equiprobable residue selection
  vector1<double> probs;
  for (unsigned i = 1; i <= pose.total_residue() - 1; ++i)  // disregard virtual residue
    probs.push_back(1);

  // Zero-out probabilities of residues that would allow folding across the cut
  invalidate_residues_spanning_cuts(pose.fold_tree(), fragments->max_frag_length(), &probs);
  numeric::normalize(probs.begin(), probs.end());

  return new BiasedFragmentMover(fragments, policy, probs);
}

void MedalAbinitioMover::estimate_residue_positions(core::pose::Pose* pose) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::fragment::FragSetOP;
  using core::scoring::ScoreFunctionFactory;
  using protocols::moves::RationalMonteCarlo;
  using protocols::nonlocal::PolicyFactory;
  assert(pose);

  FragSetOP fragments = fragments_lg_;
  RationalMonteCarlo mover(fragment_mover(*pose, fragments, PolicyFactory::get_policy("uniform", fragments)),
                           ScoreFunctionFactory::create_score_function("score0"),
                           option[OptionKeys::rigid::fragment_cycles](),
                           option[OptionKeys::abinitio::temperature](),
                           false);

  mover.apply(*pose);
}

void MedalAbinitioMover::apply(core::pose::Pose& pose) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::fragment::FragSetOP;
  using core::scoring::ScoreFunctionOP;
  using protocols::loops::Loops;
  using protocols::moves::MoverOP;
  using protocols::moves::RationalMonteCarlo;
  using protocols::moves::RigidBodyMotionMover;
  using protocols::nonlocal::PolicyOP;
  using protocols::nonlocal::PolicyFactory;
  using protocols::nonlocal::StarTreeBuilder;

  // Required options
  if (!option[OptionKeys::nonlocal::chunks].user())
    utility_exit_with_message("Failed to provide required argument -nonlocal:chunks");

  // Decompose the structure into chunks, connecting each via jump to a virtual residue
  // placed at the center of mass. The location of the virtual residue will be refined
  // as the simulation progresses.
  Loops chunks;
  decompose_structure(pose, &chunks);

  StarTreeBuilder builder;
  builder.set_up(chunks, &pose);

  // Derive reasonable estimates of residue positions through fragment insertion with vdw
  core::util::switch_to_residue_type_set(pose, core::chemical::CENTROID);
  estimate_residue_positions(&pose);

  // Set up the pose
  protocols::nonlocal::add_cutpoint_variants(&pose);
  core::scoring::constraints::add_constraints_from_cmdline_to_pose(pose);

  // Folding
  ScoreFunctionOP score = score_function();
  const double cb_start = score->get_weight(core::scoring::linear_chainbreak);
  const double cb_stop = cb_start * 2;

  Jumps jumps;
  core::pose::jumps_from_pose(pose, &jumps);
  MoverOP rigid_mover = new RigidBodyMotionMover(jumps);

  unsigned num_stages = option[OptionKeys::rigid::stages]();
  for (unsigned stage = 1; stage <= num_stages; ++stage) {
    update_position_vres(chunks, &builder, &pose);

    score->set_weight(core::scoring::linear_chainbreak,
                      numeric::linear_interpolate(cb_start, cb_stop, stage, num_stages));

    RationalMonteCarlo rigid_mc(rigid_mover, score,
                                option[OptionKeys::rigid::fragment_cycles](),
                                option[OptionKeys::abinitio::temperature](), true);

    // Large fragments in stage 1, small otherwise
    FragSetOP fragments = (stage == 1) ? fragments_lg_ : fragments_sm_;

    // Uniform policy for all stages but last
    PolicyOP policy = (stage < num_stages)
        ? PolicyFactory::get_policy("uniform", fragments)
        : PolicyFactory::get_policy("smooth", fragments);

    RationalMonteCarlo fragment_mc(fragment_mover(pose, fragments, policy), score,
                                   option[OptionKeys::rigid::fragment_cycles](),
                                   option[OptionKeys::abinitio::temperature](), true);

    // TODO(tex) consider alternative jump minimization strategies
    rigid_mc.apply(pose);
    fragment_mc.apply(pose);
  }

  // Restore simple kinematics, close remaining loops
  builder.tear_down(&pose);
  do_loop_closure(&pose);
  protocols::nonlocal::remove_cutpoint_variants(&pose);

  // Return to centroid representation and rescore
  // TODO(tex) score the pose with whatever score function you want in the output
  core::util::switch_to_residue_type_set(pose, core::chemical::CENTROID);
}

void MedalAbinitioMover::update_position_vres(const protocols::loops::Loops& chunks,
                                              protocols::nonlocal::StarTreeBuilder* builder,
                                              core::pose::Pose* pose) const {
  assert(pose);
  assert(builder);
  builder->tear_down(pose);
  builder->set_up(chunks, pose);
}

std::string MedalAbinitioMover::get_name() const {
  return "MedalAbinitioMover";
}

protocols::moves::MoverOP MedalAbinitioMover::clone() const {
  return new MedalAbinitioMover(*this);
}

protocols::moves::MoverOP MedalAbinitioMover::fresh_instance() const {
  return new MedalAbinitioMover();
}

}  // namespace medal
}  // namespace protocols
