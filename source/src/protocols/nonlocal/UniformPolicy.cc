// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/UniformPolicy.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <protocols/nonlocal/UniformPolicy.hh>

// C/C++ headers
#include <utility/assert.hh>

// Utility headers
#include <numeric/random/random.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace nonlocal {


UniformPolicy::UniformPolicy(core::fragment::FragSetCOP fragments)
    : Policy(fragments) {}

core::Size UniformPolicy::choose(const core::fragment::Frame& frame,
                                 const core::pose::Pose&) {
  assert(frame.nr_frags() > 0);
  return numeric::random::rg().random_range(1, frame.nr_frags());
}

}  // namespace nonlocal
}  // namespace protocols
