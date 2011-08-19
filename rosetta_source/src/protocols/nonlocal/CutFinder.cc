// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/CutFinder.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/CutFinder.hh>

// Utility headers
#include <numeric/random/random.hh>
#include <numeric/random/WeightedReservoirSampler.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/SecondaryStructure.hh>

// Package headers
#include <protocols/nonlocal/Region.hh>

namespace protocols {
namespace nonlocal {

static numeric::random::RandomGenerator RG(1788265078);

core::Size CutFinder::choose_cutpoint(core::Size start,
                                      core::Size stop,
                                      core::fragment::SecondaryStructureCOP ss) {
  // If secondary structure information is not available, choose a cutpoint
  // uniformly on the range [start, stop] inclusive
  if (!ss)
    return RG.random_range(start, stop);

  // Search [itvl_start, itvl_stop] for runs of consecutive loops
  utility::vector1<Region> runs;
  runs_in_range(start, stop, *ss, &runs);

  // Randomly select one of the runs in a biased manner according its length^2
  numeric::random::WeightedReservoirSampler<Region> sampler(1);
  for (core::Size i = 1; i <= runs.size(); ++i) {
    core::Size length = runs[i].length();
    sampler.consider_sample(runs[i], length * length);
  }

  // Return the midpoint of the selected run
  utility::vector1<Region> results;
  sampler.samples(&results);
  const Region& r = results[1];
  return (r.start() + r.stop()) / 2;
}

void CutFinder::runs_in_range(core::Size start,
                              core::Size stop,
                              const core::fragment::SecondaryStructure& secondary_struct,
                              utility::vector1<Region>* runs) {
  using core::Size;
  assert(runs);
  assert(start <= stop);

  int run_start = -1;
  for (Size i = start; i <= stop; ++i) {
    char ss = secondary_struct.secstruct(i);

    if (ss != 'L') {
      if (run_start >= 0) {
        runs->push_back(Region(run_start, i-1));
      }
      run_start = -1;
    } else {
      if (run_start < 0) {  // start a new run
        run_start = i;
      }
    }
  }

  // Runs are closed and their contents written to the output vector when a
  // non-loop character in the secondary structure is encountered. It's possible
  // that the interval we were given ended on a loop, in which case we need to
  // close the run and write the entry.
  if (run_start >= 0) {
    runs->push_back(Region(run_start, stop));
  }

  // In the event that the given interval contains no loop characters, randomly
  // select a position. Secondary structure prediction is not perfect, and we're
  // forced to put a cut somewhere.
  if (!runs->size()) {
    Size pos = RG.random_range(start, stop);
    runs->push_back(Region(pos, pos));
  }
}

}  // namespace nonlocal
}  // namespace protocols
