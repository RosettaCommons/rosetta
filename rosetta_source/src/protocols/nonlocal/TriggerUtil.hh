// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/TriggerUtil.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_TRIGGERUTIL_HH_
#define PROTOCOLS_NONLOCAL_TRIGGERUTIL_HH_

// External headers
#include <boost/utility.hpp>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>

namespace protocols {
namespace nonlocal {

/// @class A collection of commonly used trigger methods
class TriggerUtil : public boost::noncopyable {
 public:
  /// @brief Ramps the weight of the linear chainbreak energy term as a function
  /// of progress
  static bool ramp_chainbreaks(core::Size stage,
                               core::Size num_stages,
                               core::Size cycle,
                               core::Size num_cycles,
                               const core::pose::Pose& pose,
                               core::scoring::ScoreFunctionOP scoring);

  /// @brief Ramps the weight of the coordinate constraint energy term as a
  /// function of progress
  static bool ramp_constraints(core::Size stage,
                               core::Size num_stages,
                               core::Size cycle,
                               core::Size num_cycles,
                               const core::pose::Pose& pose,
                               core::scoring::ScoreFunctionOP scoring);

  /// @brief Ramps sequence separation as a function of progress
  static bool ramp_sequencesep(core::Size stage,
                               core::Size num_stages,
                               core::Size cycle,
                               core::Size num_cycles,
                               const core::pose::Pose& pose,
                               core::scoring::ScoreFunctionOP scoring);

  /// @brief Writes the pose to disk every N cycles
  static bool write_pose(core::Size stage,
                         core::Size num_stages,
                         core::Size cycle,
                         core::Size num_cycles,
                         const core::pose::Pose& pose,
                         core::scoring::ScoreFunctionOP scoring);

 private:
  /// @brief Computes progress. Returns a value on [0,1] inclusive.
  static double progress(core::Size p, core::Size np);

  /// @brief Updates the weight on the named scoring option to <setting>
  static void update_scoring_option(const core::scoring::ScoreType& option,
                                    core::scoring::ScoreFunctionOP scoring,
                                    core::Real setting);
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_TRIGGERUTIL_HH_
