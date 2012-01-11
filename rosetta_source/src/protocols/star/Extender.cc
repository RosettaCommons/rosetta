// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/star/Extender.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/star/Extender.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <numeric/random/random.hh>
#include <numeric/random/WeightedReservoirSampler.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/conformation/util.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/id/SequenceMapping.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <protocols/comparative_modeling/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>

namespace protocols {
namespace star {

static const double EXT_PHI = -150;
static const double EXT_PSI = +150;
static const double EXT_OMG = +180;

Extender::Extender(core::sequence::SequenceAlignmentCOP alignment, int num_residues)
    : alignment_(alignment), num_residues_(num_residues) {
  using core::id::SequenceMapping;
  using protocols::loops::Loops;
  using utility::vector1;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  assert(alignment);
  assert(num_residues > 0);

  SequenceMapping mapping = alignment->sequence_mapping(1, 2);
  vector1<int> unaligned_res;

  for (int i = 1; i < num_residues; ++i) {
    // unaligned: current residue doesn't map to any residue in the template
    // broke: current residue and next residue are both aligned, but not to consecutive residues in the template
    bool curr_aligned = mapping[i];
    bool next_aligned = mapping[i + 1];
    bool broke = (mapping[i] + 1 != mapping[i + 1]) && next_aligned;

    if (!curr_aligned) {
      unaligned_res.push_back(i);
    } else if (broke) {
      unaligned_res.push_back(i);
      unaligned_res.push_back(i+1);
    }
  }

  if (!mapping[num_residues]) {  // last residue
    unaligned_res.push_back(num_residues);
  }

  std::sort(unaligned_res.begin(), unaligned_res.end());
  utility::vector1<int>::const_iterator i = std::unique(unaligned_res.begin(), unaligned_res.end());
  unaligned_res.resize(i - unaligned_res.begin());

  int min_len = option[OptionKeys::abinitio::star::min_unaligned_len]();
  unaligned_ = protocols::comparative_modeling::pick_loops_unaligned(num_residues, unaligned_res, min_len);
  unaligned_->sequential_order();

  aligned_ = new Loops(unaligned_->invert(num_residues));
  aligned_->sequential_order();
}

void Extender::extend_section(int start, int stop, core::pose::Pose* pose) const {
  using protocols::loops::Loop;
  using protocols::loops::Loops;

  assert(pose);
  assert(start > 0);
  assert(start <= stop);
  assert(stop <= num_residues_);

  Loops region;
  region.add_loop(Loop(start, stop));
  protocols::loops::set_extended_torsions_and_idealize_loops(*pose, region);
}

void Extender::extend_unaligned(core::pose::Pose* pose) {
  using core::kinematics::FoldTree;
  using protocols::loops::Loop;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  assert(pose);
  assert(aligned_->num_loop() >= 2);

  // Keep track of the cutpoints we introduce
  cutpoints_.clear();

  const Loop& first = (*aligned_)[1];
  const Loop& last = (*aligned_)[aligned_->num_loop()];

  bool has_preceding_region = first.start() > 1;
  bool has_trailing_region = last.stop() < num_residues_;

  // beginning
  if (has_preceding_region) {
    extend_section(1, first.start() - 1, pose);
    pose->set_phi(first.start(), EXT_PHI);
    core::conformation::idealize_position(first.start(), pose->conformation());
  }

  // end
  if (has_trailing_region) {
    simple_fold_tree(pose);
    extend_section(last.stop() + 1, num_residues_, pose);
    pose->set_omega(last.stop(), EXT_OMG);
    core::conformation::idealize_position(last.stop(), pose->conformation());
  }

  // interior
  for (unsigned i = 2; i <= aligned_->num_loop(); ++i) {
    const Loop& prev = (*aligned_)[i - 1];
    const Loop& curr = (*aligned_)[i];

    // If the sequence separation between consecutive aligned regions <= x,
    // chain the unaligned residues off the first aligned region. Otherwise,
    // split the unaligned residues between the two aligned regions.
    unsigned separation = curr.start() - prev.stop();

    if (separation <= option[OptionKeys::abinitio::star::short_loop_len]()) {
      FoldTree tree(num_residues_);
      tree.new_jump(1, curr.start(), curr.start() - 1);
      pose->fold_tree(tree);
      extend_section(prev.stop(), curr.start() - 2, pose);

      // Since we chained the unaligned residues off one side (i.e. `prev`),
      // set the cutpoint before the first residue on the other side (i.e. `curr`)
      cutpoints_.push_back(curr.start() - 1);
    } else {
      unsigned cut = choose_cutpoint(prev.stop() + 1, curr.start() - 1);
      cutpoints_.push_back(cut);

      FoldTree tree(num_residues_);
      tree.new_jump(1, curr.start() + 1, cut);
      pose->fold_tree(tree);

      // Extend torsions left of the cut-- incr order toward cut
      for (unsigned i = prev.stop(); i < cut; ++i) {
        pose->set_phi(i, EXT_PHI);
        pose->set_psi(i, EXT_PSI);
        pose->set_omega(i, EXT_OMG);
        core::conformation::idealize_position(i, pose->conformation());
      }

      // Extend torsions right of the cut-- decr order toward cut
      for (unsigned i = curr.start(); i >= cut; --i) {
        pose->set_phi(i, EXT_PHI);
        pose->set_psi(i, EXT_PSI);
        pose->set_omega(i, EXT_OMG);
        core::conformation::idealize_position(i, pose->conformation());
      }
    }
  }

  simple_fold_tree(pose);
}

int Extender::choose_cutpoint(int start, int stop) const {
  using numeric::random::WeightedReservoirSampler;
  using utility::vector1;

  assert(start > 0);
  assert(start <= stop);

  if (!pred_ss_) {
    return numeric::random::random_range(start, stop);
  }

  WeightedReservoirSampler<int> sampler(1);

  for (int i = start; i <= stop; ++i) {
    double weight = pred_ss_->loop_fraction(i) + 1e-20;
    sampler.consider_sample(i, weight);
  }

  vector1<int> samples;
  sampler.samples(&samples);

  assert(sampler.num_considered());
  return samples[1];
}

void Extender::simple_fold_tree(core::pose::Pose* pose) const {
  assert(pose);
  core::kinematics::FoldTree tree(num_residues_);
  pose->fold_tree(tree);
}

}  // namespace star
}  // namespace protocols
