// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/CutFinder.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_CUTFINDER_HH_
#define PROTOCOLS_NONLOCAL_CUTFINDER_HH_

// Package headers
#include <protocols/nonlocal/Region.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace nonlocal {

/// @class Contains utility methods for identifying preferable cutpoint
/// locations within non-local groupings. Biases selection toward runs of
/// consecutive loop characters, as predicted by secondary structure.
class CutFinder {
 public:
  /// @brief Randomly select cutpoints between adjacent non-local fragments.
	/// If <ss> is not null, use it to improve cutpoint selection.
  static core::Size choose_cutpoint(core::Size start,
                                    core::Size stop,
                                    core::fragment::SecondaryStructureCOP ss);

 private:
  /// @brief Searches <predicted_secondary_struct_> for runs of consecutive
  /// loop characters (L). Tracks the starting position and length of each run.
  /// Stores the result in the output parameter <runs>.
  static void runs_in_range(core::Size start,
                            core::Size stop,
                            const core::fragment::SecondaryStructure& ss,
                            utility::vector1<Region>* runs);
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_CUTFINDER_HH_
