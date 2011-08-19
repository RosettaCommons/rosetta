// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/TreeBuilder.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_TREEBUILDER_HH_
#define PROTOCOLS_NONLOCAL_TREEBUILDER_HH_

// Unit headers
#include <protocols/nonlocal/TreeBuilder.fwd.hh>

// Package headers
#include <protocols/nonlocal/NLGrouping.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/fragment/SecondaryStructure.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace nonlocal {

class TreeBuilder : public utility::pointer::ReferenceCount {
 public:
  /// @brief Modifies <pose>'s FoldTree given non-local contact information.
  /// Secondary structure, if provided, is used to improve cutpoint selection.
  virtual void build(const NLGrouping& grouping,
                     core::pose::Pose* pose,
                     core::fragment::SecondaryStructureOP secondary_struct = NULL) = 0;

  /// @brief Reverts any modifications to <pose> that this builder may have
  // introduced in previous calls to build(). Only subclasses that introduce
  // modifications are responsible for overriding this method.
  virtual void tear_down(core::pose::Pose* pose) {}
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_TREEBUILDER_HH_
