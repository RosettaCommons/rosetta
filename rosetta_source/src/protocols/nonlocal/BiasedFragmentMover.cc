// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/BiasedFragmentMover.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/nonlocal/BiasedFragmentMover.hh>

// C/C++ headers
#include <algorithm>
#include <string>

// External headers
#include <boost/format.hpp>

// Utility headers
#include <basic/Tracer.hh>
#include <numeric/prob_util.hh>
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>

// Package headers
#include <protocols/nonlocal/Policy.hh>

namespace protocols {
namespace nonlocal {

typedef core::fragment::FragSetOP FragSetOP;
typedef utility::vector1<double> Probabilities;

static basic::Tracer TR("protocols.nonlocal.BiasedFragmentMover");

BiasedFragmentMover::BiasedFragmentMover(const FragSetOP& fragments, const PolicyOP& policy, const Probabilities& pdf)
    : fragments_(fragments), policy_(policy), pdf_(pdf) {
  assert(fragments);
  assert(policy);
  initialize_library();
  initialize_probabilities();
}

/// @detail Creates a position-indexable Frame lookup
void BiasedFragmentMover::initialize_library() {
  unsigned position;
  for (core::fragment::FrameIterator i = fragments_->begin(); i != fragments_->end(); ++i) {
    position = (*i)->start();
    frames_[position] = **i;
  }
}

/// @detail Normalizes the pdf, then computes the cdf
void BiasedFragmentMover::initialize_probabilities() {
  numeric::normalize(pdf_.begin(), pdf_.end());
  std::copy(pdf_.begin(), pdf_.end(), cdf_.begin());
  numeric::cumulative(cdf_.begin(), cdf_.end());
}

/// @detail Verifies that the probability of selecting invalid positions is 0
void BiasedFragmentMover::verify_probabilities_or_die(const core::kinematics::FoldTree& tree) const {
  const unsigned fragment_len = fragments_->max_frag_length();

  for (unsigned i = 1; i <= tree.num_cutpoint(); ++i) {
    const unsigned cutpoint = tree.cutpoint(i);

    // In order to avoid folding across the cut, certain residues must have zero probability of selection
    for (unsigned j = (cutpoint - fragment_len + 2); j <= cutpoint; ++j) {
      const double prob = pdf_[j];
      if (prob > 0) {
        utility_exit_with_message((boost::format("Non-zero selection probability for unmodifiable residue %1%") % j).str());
      }
    }
  }
}

void BiasedFragmentMover::apply(core::pose::Pose& pose) {
  using core::fragment::Frame;
  using core::kinematics::MoveMap;

  // Verify user-specified sampling probabilities
  verify_probabilities_or_die(pose.fold_tree());

  const unsigned insertion_pos = random_position();
  const Frame& frame = frames_[insertion_pos];
  const unsigned fragment_num = policy_->choose(frame, pose);

  MoveMap movable;
  movable.set_bb(true);
  frame.apply(movable, fragment_num, pose);
  TR.Debug << "Inserted fragment at position " << insertion_pos << std::endl;
}

/// @detail Selects the insertion position in a weighted random fashion using binary search on the cdf
unsigned BiasedFragmentMover::random_position() const {
  Probabilities::const_iterator i = std::lower_bound(cdf_.begin(), cdf_.end(), numeric::random::uniform());
  return i - cdf_.begin() + 1;
}

std::string BiasedFragmentMover::get_name() const {
  return "BiasedFragmentMover";
}

}  // namespace nonlocal
}  // namespace protocols
