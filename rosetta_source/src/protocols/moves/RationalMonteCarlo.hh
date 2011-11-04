// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/RationalMonteCarlo.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_MOVES_RATIONAL_MONTE_CARLO_HH_
#define PROTOCOLS_MOVES_RATIONAL_MONTE_CARLO_HH_

// Unit header
#include <protocols/moves/RationalMonteCarlo.fwd.hh>

// C/C++ headers
#include <string>

// External headers
#include <boost/function.hpp>
#include <boost/unordered/unordered_map.hpp>

// Project headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Package headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/moves/MonteCarlo.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace moves {

/// @brief Trigger API definition
typedef boost::function<void(const core::pose::Pose&)> RationalMonteCarloTrigger;
typedef boost::unordered_map<int, RationalMonteCarloTrigger> Triggers;

/// @class Trial-based Monte Carlo minization primitive. Do not modify this
/// class; almost anything you could possibly want to add is a bad idea.
class RationalMonteCarlo : public protocols::moves::Mover {
  typedef core::Real Real;
  typedef core::Size Size;
  typedef core::pose::Pose Pose;
  typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;

 public:
  RationalMonteCarlo(MoverOP mover,
                     ScoreFunctionOP score,
                     Size num_trials,
                     Real temperature,
                     bool recover_low);

  /// @brief Applies the underlying mover to <pose> the specified number of times
  void apply(Pose& pose);

  /// @brief Returns the number of trials to be attempted in calls to apply()
  Size num_trials() const;

  /// @brief Returns true if we should recover the low scoring pose at the end
  /// of apply, false otherwise.
  bool recover_low() const;

  /// @brief Returns this mover's name
  std::string get_name() const;

  /// @brief Registers the specified trigger with this instance. Returns a unique
  /// identifier for referring to this trigger for subsequent operations (e.g. remove).
  int add_trigger(const RationalMonteCarloTrigger& trigger);

  /// @brief Unregisters the specified trigger.
  void remove_trigger(int trigger_id);

 protected:
  /// @brief Executes all triggers attached to this instance. The order of trigger
  /// execution is undefined. Do not assume, depend, or in any way rely upon a
  /// partiular ordering.
  void fire_all_triggers(const Pose& pose);

 private:
  /// @brief Underlying mover
  MoverOP mover_;

  /// @brief Determines whether applications of the base mover should be accepted
  /// by applying the Metropolis criterion to the pose
  MonteCarloOP mc_;

  /// @brief Number of times to execute the underlying mover in calls to apply()
  Size num_trials_;

  /// @brief Determines whether the low scoring pose should be recovered at the
  /// conclusion of the apply() method
  bool recover_low_;

  /// @brief Collection of function callbacks
  Triggers triggers_;

  /// @brief Next trigger id to be assigned
  Size next_trigger_id_;
};

}  // namespace moves
}  // namespace protocols

#endif  // PROTOCOLS_MOVES_RATIONAL_MONTE_CARLO_HH_
