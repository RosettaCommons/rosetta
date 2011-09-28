// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/medal/MedalMover.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_MEDAL_MEDAL_MOVER_HH_
#define PROTOCOLS_MEDAL_MEDAL_MOVER_HH_

// Unit header
#include <protocols/medal/MedalMover.fwd.hh>

// C/C++ headers
#include <string>

// External headers
#include <boost/unordered/unordered_map.hpp>

// Project headers
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace medal {

class MedalMover : public protocols::moves::Mover {
  typedef boost::unordered_map<int, core::kinematics::Jump> Jumps;

 public:
  void apply(core::pose::Pose& pose);

  // -- jd2 -- //
  std::string get_name() const;
  protocols::moves::MoverOP clone() const;
  protocols::moves::MoverOP fresh_instance() const;

private:
  /// @brief Closes chainbreaks in <pose>
  void do_loop_closure(core::pose::Pose* pose) const;

  /// @brief Performs kinematically-aware, scored fragment insertion
  void do_fragment_insertion(const core::scoring::ScoreFunctionOP& score,
                             core::pose::Pose* pose) const;

  /// @brief Performs rigid body moves
  void do_rigid_body_moves(const core::scoring::ScoreFunctionOP& score,
                           core::pose::Pose* pose) const;

  /// @brief Retrieves jump information from <pose>, storing the result in <jumps>
  void jumps_from_pose(const core::pose::Pose& pose, Jumps* jumps) const;

  /// @brief Configures a basic score functions which callers can then specialize
  core::scoring::ScoreFunctionOP base_score_function() const;

  /// @brief Configures the score function used for rigid body perturbation
  core::scoring::ScoreFunctionOP perturb_score_function() const;

  /// @brief Configures the score function used for fragment insertion
  core::scoring::ScoreFunctionOP insert_score_function() const;
};

}  // namespace medal
}  // namespace protocols

#endif  // PROTOCOLS_MEDAL_MEDAL_MOVER_HH_
