// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/BoundaryFinder.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_BOUNDARYFINDER_HH_
#define PROTOCOLS_NONLOCAL_BOUNDARYFINDER_HH_

// Project headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>

namespace protocols {
namespace nonlocal {

/// @class Contains utility methods for identifying the cut- and jump-point
/// boundaries on both sides of a residue by examining the fold tree.
class BoundaryFinder {
  typedef core::Size Size;

 public:
  /// @brief Given a valid residue position, searches the fold tree to identify
  /// the most immediate upper and lower boundaries. Stores the result in output
  /// parameters <upper> and <lower>.
  static void boundaries(const core::kinematics::FoldTree& tree,
                         Size position,
                         Size* lower,
                         Size* upper);
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_BOUNDARYFINDER_HH_
