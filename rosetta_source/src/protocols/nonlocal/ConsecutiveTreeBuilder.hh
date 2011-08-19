// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/ConsecutiveTreeBuilder.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_CONSECUTIVE_TREE_BUILDER_HH_
#define PROTOCOLS_NONLOCAL_CONSECUTIVE_TREE_BUILDER_HH_

// Unit headers
#include <protocols/nonlocal/ConsecutiveTreeBuilder.fwd.hh>

// Package headers
#include <protocols/nonlocal/NLGrouping.hh>
#include <protocols/nonlocal/TreeBuilder.hh>

// Project headers
#include <core/fragment/SecondaryStructure.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace nonlocal {

class ConsecutiveTreeBuilder : public TreeBuilder {
 public:
  void build(const NLGrouping& groupings,
						 core::pose::Pose* pose,
             core::fragment::SecondaryStructureOP secondary_struct = 0);
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_CONSECUTIVE_TREE_BUILDER_HH_
