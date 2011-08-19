// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/BoundaryFinder.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/nonlocal/BoundaryFinder.hh>

// C/C++ headers
#include <cassert>

// Project headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>

namespace protocols {
namespace nonlocal {

void BoundaryFinder::boundaries(const core::kinematics::FoldTree& tree,
                                core::Size position,
                                core::Size* lower,
                                core::Size* upper) {
  using core::Size;

  assert(lower);
  assert(upper);

  Size lower_cut = 0;
  Size upper_cut = tree.nres();
  Size num_cutpoints = tree.num_cutpoint();
  for (Size i = 1; i <= num_cutpoints; ++i) {
    Size cutpoint = tree.cutpoint(i);

    // find the upper boundary (inclusive)
    if (cutpoint >= position && cutpoint < upper_cut)
      upper_cut = cutpoint;

    // find the lower boundary (exclusive)
    if (cutpoint < position && cutpoint > lower_cut)
      lower_cut = cutpoint;
  }

  // set output parameters
  *lower = lower_cut + 1;
  *upper = upper_cut;
}

}  // namespace nonlocal
}  // namespace protocols
