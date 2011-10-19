// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/RigidBodyMotionMover.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_MOVES_RIGID_BODY_MOTION_MOVER_HH_
#define PROTOCOLS_MOVES_RIGID_BODY_MOTION_MOVER_HH_

// Unit header
#include <protocols/moves/RigidBodyMotionMover.fwd.hh>

// C/C++ headers
#include <set>
#include <string>

// Project headers
#include <core/kinematics/Jump.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace moves {

class RigidBodyMotionMover : public Mover {
  typedef std::set<int> Jumps;

 public:
  explicit RigidBodyMotionMover(const Jumps& jumps);

  /// @brief Randomly selects a chunk and perturbs its rigid body transformation
  void apply(core::pose::Pose& pose);

  // jd2
  std::string get_name() const;
  MoverOP clone() const;

 private:
	/// @brief Returns a randomly selected jump
	int random_jump() const;

  Jumps jumps_;
};

}  // namespace moves
}  // namespace protocols

#endif  // PROTOCOLS_MOVES_RIGID_BODY_MOTION_MOVER_HH_
