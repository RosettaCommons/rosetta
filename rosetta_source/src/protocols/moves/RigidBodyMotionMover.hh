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
#include <core/pose/Pose.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace moves {

class RigidBodyMotionMover : public Mover {
  typedef std::set<int> Jumps;

 protected:
  static const int TRANS_X = 1;
  static const int TRANS_Y = 2;
  static const int TRANS_Z = 3;
  static const int ROT_X = 4;
  static const int ROT_Y = 5;
  static const int ROT_Z = 6;

 public:
  RigidBodyMotionMover();
  explicit RigidBodyMotionMover(const Jumps& jumps);

  /// @brief Randomly selects a chunk and perturbs its rigid body transformation
  void apply(core::pose::Pose& pose);

  /// @brief Returns the magnitude of translation
  double magnitude_translation() const;

  /// @brief Returns the magnitude of rotation
  double magnitude_rotation() const;

  /// @brief Updates the magnitude of translation
  void set_magnitude_translation(double mag_trans);

  /// @brief Updates the magnitude of rotation
  void set_magnitude_rotation(double mag_rot);

  /// @brief Returns the name of this mover
  std::string get_name() const;

  /// @brief Returns a new instance created by the copy constructor
  MoverOP clone() const;

  /// @brief Returns a new instance created by the default constructor
  MoverOP fresh_instance() const;

 private:
  /// @brief Shared initialization
  void initialize(const Jumps& jumps);

  /// @brief Returns a randomly selected jump
  int random_jump() const;

  /// @brief Collection of jumps to manipulate
  Jumps jumps_;

  /// @brief Magnitude of rotation. Unless otherwise specified, equal to -rigid:rotation.
  double mag_rot_;

  /// @brief Magnitude of translation. Unless otherwise specified, equal to -rigid:translation.
  double mag_trans_;
};

}  // namespace moves
}  // namespace protocols

#endif  // PROTOCOLS_MOVES_RIGID_BODY_MOTION_MOVER_HH_
