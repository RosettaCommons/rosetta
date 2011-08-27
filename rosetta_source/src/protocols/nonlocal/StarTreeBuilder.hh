// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/StarTreeBuilder.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_STARTREEBUILDER_HH_
#define PROTOCOLS_NONLOCAL_STARTREEBUILDER_HH_

// Unit headers
#include <protocols/nonlocal/StarTreeBuilder.fwd.hh>

// Package headers
#include <protocols/nonlocal/NLGrouping.hh>
#include <protocols/nonlocal/TreeBuilder.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace nonlocal {

class StarTreeBuilder : public TreeBuilder {
 public:
  /// @brief Assigns <virtual_res_> a value signifying that it is uninitialized
  StarTreeBuilder();

	/// @brief Constructs a star fold tree by placing a virtual residue at
	/// <grouping>'s center of mass and adding jumps to a stochastically
	/// chosen residue in each chunk.
  void set_up(const NLGrouping& grouping, core::pose::Pose* pose);

  /// @brief Removes the virtual residue added to <pose> in calls to build()
  void tear_down(core::pose::Pose* pose);

 private:
  /// @brief Stochastically selects an anchor position on [fragment.start(), fragment.stop()]
  /// according to per-residue structural conservation
  core::Size choose_conserved_position(const NLFragmentGroup& fragment,
                                       const core::pose::Pose& pose);

  /// @brief Returns the index of the virtual residue placed at the center of
  /// mass in build(). If no such call has occurred, returns -1.
  int virtual_residue() const;

  /// @brief Sets the index of the virtual residue
  void virtual_residue(int virtual_res);

  /// @brief Index of the virtual residue we added to the pose in build()
  int virtual_res_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_STARTREEBUILDER_HH_
