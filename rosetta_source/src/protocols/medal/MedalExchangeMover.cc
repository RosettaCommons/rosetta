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
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <numeric/prob_util.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/util/kinematics_util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/comparative_modeling/ThreadingJob.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/medal/util.hh>
#include <protocols/nonlocal/BiasedFragmentMover.hh>
#include <protocols/nonlocal/Policy.hh>
#include <protocols/nonlocal/PolicyFactory.hh>
#include <protocols/nonlocal/TreeBuilder.hh>
#include <protocols/nonlocal/TreeBuilderFactory.hh>
#include <protocols/nonlocal/util.hh>
#include <protocols/simple_moves/rational_mc/RationalMonteCarlo.hh>
#include <protocols/star/Extender.hh>

namespace protocols {
namespace medal {

using core::Real;
using core::Size;
using core::scoring::constraints::ConstraintSet;
using core::scoring::constraints::ConstraintSetOP;
using core::scoring::constraints::ConstraintSetCOP;
using protocols::loops::Loop;
using protocols::loops::Loops;
using protocols::loops::LoopsOP;
using protocols::loops::LoopsCOP;
using utility::vector1;

static basic::Tracer TR("protocols.medal.MedalExchangeMover");

void to_centroid(core::pose::Pose* pose) {
  assert(pose);
  if (!pose->is_centroid())
    core::util::switch_to_residue_type_set(*pose, core::chemical::CENTROID);
}

/// @detail Combines both sets of loops, sorting the result in increasing order of start position
protocols::loops::LoopsCOP combine_loops(LoopsCOP aligned, LoopsCOP unaligned) {
  assert(aligned);
  assert(unaligned);

  LoopsOP combined = new Loops();

  for (Loops::const_iterator i = aligned->begin(); i != aligned->end(); ++i)
    combined->add_loop(*i);

  for (Loops::const_iterator i = unaligned->begin(); i != unaligned->end(); ++i)
    combined->add_loop(*i);

  combined->sequential_order();
  return combined;
}

void MedalExchangeMover::setup_constraints(const core::pose::Pose& pose, LoopsCOP aligned, ConstraintSetOP constraints) const {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::PointPosition;
  using core::conformation::Residue;
  using core::id::AtomID;
  using core::scoring::constraints::BoundFunc;
  using core::scoring::constraints::ConstraintCOP;
  using core::scoring::constraints::CoordinateConstraint;
  using core::scoring::constraints::FuncOP;
  assert(aligned);
  assert(constraints);

  const Real delta = option[OptionKeys::cm::sanitize::bound_delta]();
  const Real sd = option[OptionKeys::cm::sanitize::bound_sd]();

  // Fixed reference position
  const AtomID fixed_atom(1, pose.total_residue());
  const PointPosition& fixed_coords = pose.xyz(fixed_atom);

  for (Size i = 1; i <= aligned->size(); ++i) {
    const Loop& region = (*aligned)[i];

    for (Size j = region.start(); j <= region.stop(); ++j) {
      const Residue& residue = pose.conformation().residue(j);
      const AtomID ca_atom(residue.atom_index("CA"), j);
      const PointPosition& ca_coords = pose.xyz(ca_atom);

      Real distance = ca_coords.distance(fixed_coords);
      FuncOP function = new BoundFunc(distance - delta, distance + delta, sd, "");
      ConstraintCOP constraint = new CoordinateConstraint(ca_atom, fixed_atom, fixed_coords, function);
      constraints->add_constraint(constraint);
    }
  }

  constraints->show(TR.Debug);
  TR.flush_all_channels();
}

void MedalExchangeMover::setup_sampling_probs(Size num_residues, const core::kinematics::FoldTree& tree, LoopsCOP aligned, vector1<double>* probs) const {
  assert(probs);
  assert(aligned);

  probs->resize(num_residues, 0);

  for (Loops::const_iterator i = aligned->begin(); i != aligned->end(); ++i) {
    for (Size j = i->start(); j <= i->stop(); ++j) {
      (*probs)[j] = 1;
    }
  }

  protocols::medal::invalidate_residues_spanning_cuts(tree, fragments_->max_frag_length(), probs);
  numeric::normalize(probs->begin(), probs->end());
  numeric::print_probabilities(*probs, TR.Debug);
}

void MedalExchangeMover::apply(core::pose::Pose& pose) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::id::AtomID;
  using core::scoring::Energies;
  using core::scoring::ScoreFunction;
  using core::scoring::ScoreFunctionOP;
  using protocols::comparative_modeling::ThreadingJob;
  using protocols::moves::MoverOP;
  using protocols::nonlocal::BiasedFragmentMover;
  using protocols::nonlocal::PolicyOP;
  using protocols::nonlocal::PolicyFactory;
  using protocols::nonlocal::TreeBuilderOP;
  using protocols::nonlocal::TreeBuilderFactory;
  using protocols::simple_moves::rational_mc::RationalMonteCarlo;
  using protocols::star::Extender;

  ThreadingJob const * const job = protocols::nonlocal::current_job();
  to_centroid(&pose);

  Extender extender(job->alignment().clone(), pose.total_residue());
  extender.set_secondary_structure(pred_ss_);
  extender.extend_unaligned(&pose);

  LoopsCOP aligned = extender.aligned();
  LoopsCOP unaligned = extender.unaligned();
  LoopsCOP combined = combine_loops(aligned, unaligned);
  TR << "Aligned:" << std::endl << *aligned << std::endl;

  TreeBuilderOP builder = TreeBuilderFactory::get_builder("star");
  builder->set_up(*combined, &pose);
  TR << pose.fold_tree() << std::endl;

  ConstraintSetOP constraints = new ConstraintSet();
  setup_constraints(pose, aligned, constraints);
  pose.constraint_set(constraints);

  ScoreFunctionOP score = new ScoreFunction();
  score->set_weight(core::scoring::coordinate_constraint, 1);

  // Sampling parameters
  const Real temp = option[OptionKeys::abinitio::temperature]();
  const Size num_fragments = option[OptionKeys::cm::sanitize::num_fragments]();
  const Size num_cycles = static_cast<Size>(aligned->nr_residues() * 100 * option[OptionKeys::abinitio::increase_cycles]());

  vector1<double> probs;
  setup_sampling_probs(pose.total_residue() - 1, pose.fold_tree(), aligned, &probs);

  PolicyOP policy = PolicyFactory::get_policy("uniform", fragments_, num_fragments);
  MoverOP fragment_mover = new BiasedFragmentMover(policy, probs);
  RationalMonteCarlo mover(fragment_mover, score, num_cycles, temp, false);
  mover.apply(pose);

  // Removing the virtual residue introduced by the star fold tree invalidates
  // the pose's cached energies. Doing so causes the score line in the silent
  // file to be 0.
  Energies energies = pose.energies();

  // Housekeeping
  pose.remove_constraints();
  builder->tear_down(&pose);

  pose.set_new_energies_object(new Energies(energies));
}

MedalExchangeMover::MedalExchangeMover() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using core::fragment::FragmentIO;
  using core::fragment::FragSetOP;
  using core::fragment::SecondaryStructure;

  FragmentIO io;
  fragments_ = io.read_data(option[in::file::frag3]());
  pred_ss_ = new SecondaryStructure(*fragments_);
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
