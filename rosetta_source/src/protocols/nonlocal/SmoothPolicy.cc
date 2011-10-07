// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/SmoothPolicy.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/nonlocal/SmoothPolicy.hh>

// Utility headers
#include <numeric/random/random.hh>
#include <utility/minmax.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace nonlocal {

typedef SmoothPolicy::Candidate Candidate;
typedef utility::vector1<core::Real> ScoreList;
typedef utility::vector1<Candidate> CandidateList;

SmoothPolicy::SmoothPolicy(core::fragment::FragSetCOP fragments)
    : Policy(fragments) {}

core::Size SmoothPolicy::choose(const core::fragment::Frame& frame,
                                const core::pose::Pose& pose) {
  assert(frame.nr_frags() > 0);

  ScoreList scores;
  scorer_.score(frame, pose, scores);

  CandidateList candidates;
  for (core::Size i = 1; i <= scores.size(); ++i) {
    core::Real score = scores[i];
    if (score < scorer_.cutoff()) {
      candidates.push_back(Candidate(score, i));
    }
  }

  // If no candidates met the score threshold, return the best scoring fragment.
  // Otherwise, uniformly choose among the candidates.
  if (candidates.size() == 0) {
    return utility::argmin(scores);
  } else {
    return numeric::random::random_range(1, candidates.size());
  }
}

}  // namespace nonlocal
}  // namespace protocols
