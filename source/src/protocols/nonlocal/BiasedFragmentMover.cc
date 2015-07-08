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
#include <iterator>
#include <string>

// External headers
#include <boost/format.hpp>

// Utility headers
#include <numeric/prob_util.hh>
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

#include <core/fragment/FrameIteratorWorker_.hh>

//Auto Headers
#include <numeric/random/random.fwd.hh>

namespace protocols {
namespace nonlocal {

typedef utility::vector1<double> Probabilities;

BiasedFragmentMover::BiasedFragmentMover(const PolicyOP& policy, const Probabilities& pdf)
    : policy_(policy), pdf_(pdf) {
  assert(policy);
  fragments_ = policy->fragments();
  initialize_library();
  initialize_probabilities();

  // Avoid repeated MoveMap construction in apply()
  movable_.set_bb(true);
}

/// @detail Creates a position-indexable Frame lookup
void BiasedFragmentMover::initialize_library() {
  //unsigned position;
  for (core::fragment::ConstFrameIterator i = fragments_->begin(); i != fragments_->end(); ++i) {
    unsigned position = (*i)->start();
    frames_[position] = **i;
  }
}

/// @detail Computes cdf from pdf. cumulative() takes care of normalization.
void BiasedFragmentMover::initialize_probabilities() {
  std::copy(pdf_.begin(), pdf_.end(), std::back_inserter(cdf_));
  numeric::cumulative(cdf_.begin(), cdf_.end());
}

/// @detail Verifies that the probability of selecting invalid positions is 0
void BiasedFragmentMover::verify_probabilities_or_die(const core::kinematics::FoldTree& tree) const {
  const unsigned fragment_len = fragments_->max_frag_length();

  for (int i = 1; i <= tree.num_cutpoint(); ++i) {
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
  // Verify user-specified sampling probabilities
  verify_probabilities_or_die(pose.fold_tree());

  // Select insertion position and fragment
  const unsigned insertion_pos = random_position();
  const core::fragment::Frame& frame = frames_[insertion_pos];
  const unsigned fragment_num = policy_->choose(frame, pose);

  frame.apply(movable_, fragment_num, pose);
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
